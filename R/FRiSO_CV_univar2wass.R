#' K-fold cross-validation (CV) for FRiSO, for univariate distributions with 2-Wasserstein metric
#'
#' @description
#' This function performs K-fold cross-validation with the Fréchet Ridge
#' Selection Operator (FRiSO; \insertCite{tucker_variable_2023}{fastfrechet}),
#' for the space of univariate distribution responses equipped with the
#' 2-Wasserstein metric. The optimal hyperparameter \eqn{\tau} is chosen as that
#' which minimizes out-of-sample error among a collection of user-supplied
#' values in `tauseq`. If there are ties for the minimum error, the model formed
#' from the smallest \eqn{\tau} is chosen. The model corresponding to the
#' optimal \eqn{\tau} is also provided as part of the output list object.
#'
#' @inheritParams frechetreg_univar2wass
#' @param K Numeric whole number no bigger than `nrow(X)`.
#' @param thresh Positive scalar that determines the selection cutoff for the `lambda` vector in the variable selection problem.
#' @param ... other inputs to [fastfrechet::FRiSO_univar2wass()].
#'
#' @return A list object with components:
#' \tabular{ll}{
#'   `tauseq` \tab returns the numeric vector containing \eqn{\tau} values at which FRiSO problem was solved. \cr
#'   `errors` \tab returns a (`p` \eqn{\times} `length(tauseq)`) matrix that contains the refitted test errors.  \cr
#'   `error_sum` \tab returns a `length(tauseq)`-long vector that contains the aggregate CV error for each \eqn{\tau}. \cr
#'   `opt_tau` \tab returns a numeric scalar that is the \eqn{\tau} that minimizes the aggregate CV error. \cr
#'   `opt_lambda` \tab returns a `p`-long numeric vector that is the high-precision 'allowance vector'. \cr
#'   `opt_selected` \tab returns a numeric vector (up to `p`-long) that contains the indices of the variables selected in the optimal model. \cr
#' }
#'
#' @references
#' \insertRef{tucker_variable_2023}{fastfrechet}
#'
#' @export
#'
#' @examples
#' # Generate data for X and Y inputs by using the output of `generate_zinbinom_qf`
#' n <- 100 # number of samples - nrow(X) and nrow(Y).
#' p <- 10 # number of covariates - ncol(X).
#' m <- 50 # EQF grid density - ncol(Y).
#' lower <- 0
#' upper <- Inf
#'
#' set.seed(31)
#' mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
#' X <- mydata$X # (n x p) matrix of covariates
#' Y <- mydata$Y # (n x m) matrix of EQFs, stored row-wise
#'
#' # Set cross-validation and FRiSO parameters
#' K <- 10
#' thresh <- 0.0001
#' tauseq <- seq(0.2, 20, 0.2)
#' nudge <- 0.001
#' eps <- 0.001
#'
#' # Run cross-validation:
#' cv <- FRiSO_CV_univar2wass(
#'   X = X,
#'   Y = Y,
#'   K = K,
#'   thresh = thresh,
#'   lower = lower,
#'   upper = upper,
#'   tauseq = tauseq,
#'   nudge = nudge,
#'   eps = eps
#' )
#'
#' # Plot errors per fold and average fold error:
#' matplot(tauseq, cv$errors, type = "l", lty = 1, main = "CV Fold Errors")
#' lines(tauseq, cv$error_sum / K, lwd = 3)
#' points(cv$opt_tau, min(cv$error_sum) / K, pch = 1, lwd = 2, cex = 1.5)
#'
#' # Identify which variables are selected in "optimal" model:
#' cv$opt_selected
FRiSO_CV_univar2wass <- function(X,
                                 Y,
                                 K = NULL,
                                 thresh = 0.0001,
                                 ...) {
  # Extract call parameters:
  Call <- c(as.list(environment()), list(...))

  # Create call list for `FRiSO_univar2wass`, matching from Call if provided:
  Send_Call <- formals(FRiSO_univar2wass)
  Send_Call[names(Call)[names(Call) %in% names(Send_Call)]] <- Call[names(Call) %in% names(Send_Call)]

  # Extract necessary parameters for `frechetreg_univar2wass`:
  full_tauseq <- Send_Call$"tauseq"
  lower <- Send_Call$"lower"
  upper <- Send_Call$"upper"

  # Check numeric, and grab dimensions:
  check_numeric(X, "matrix", finite = TRUE)
  check_numeric(Y, "matrix", finite = TRUE)
  n <- nrow(X)
  m <- ncol(Y)
  p <- ncol(X)

  # Dimension checks:
  if (n != nrow(Y)) stop("X and Y must have same number of rows.")

  # If K is not specified, then we set it equal to n, i.e. LOOCV:
  if (is.null(K)) {
    K <- n
  } else {
    check_wholenumber(K)
    if (K > n) stop("'K' must be a whole number no greater than nrow(X).")
  }

  # Check numeric for thresh:
  check_numeric(thresh, "scalar", finite = TRUE)
  if (thresh <= 0) stop("'thresh' must be a positive scalar.")

  # Subset into K folds:
  K_sub <- sample(rep(1:K, ceiling(n / K))[1:n])

  # Create empty CV errors array:
  errors <- matrix(NA, nrow = length(full_tauseq), ncol = K)

  # Loop through folds:
  for (k in 1:K) {
    # Collect training and test data:
    Xtrain <- X[K_sub != k, , drop = FALSE]
    Ytrain <- Y[K_sub != k, , drop = FALSE]

    Xtest <- X[K_sub == k, , drop = FALSE]
    Ytest <- Y[K_sub == k, , drop = FALSE]

    # Train, i.e. run FRiSO:
    Send_Call$"X" <- Xtrain
    Send_Call$"Y" <- Ytrain
    L <- do.call(FRiSO_univar2wass, args = Send_Call)

    # Calculate errors from test subset:
    C_init <- matrix(0, nrow(Xtest), m + 1)
    for (t in 1:length(full_tauseq)) {
      # Find chosen model support:
      selected <- which(L[, t] > thresh)

      # If no variables are selected, i.e. 'thresh' is very large compared to
      # min('tauseq') / ncol(X), then train on all-zeros matrix:
      if (length(selected) == 0) {
        output <- frechetreg_univar2wass(
          X = matrix(0, nrow(Xtrain), 1),
          Y = Ytrain,
          Z = matrix(0, nrow(Xtest), 1),
          C_init = C_init,
          lambda = NULL,
          lower = lower,
          upper = upper
        )
      } else {
        output <- frechetreg_univar2wass(
          X = Xtrain[, selected, drop = FALSE],
          Y = Ytrain,
          Z = Xtest[, selected, drop = FALSE],
          C_init = C_init,
          lambda = NULL,
          lower = lower,
          upper = upper
        )
      }

      # Save Lagrange multiplier for next warm start:
      C_init <- output$"Lagrange_Multiplier"

      # Calculate refitted test errors:
      errors[t, k] <- sum((output$"Qhat" - Ytest)^2) / m
    }
  }

  # Calculate aggregate CV error for each \tau:
  error_sum <- rowSums(errors)

  # Calculate \tau which minimizes aggregate CV error:
  opt_tau <- full_tauseq[which.min(error_sum)[1]]

  # Calculate high-precision 'allowance vector' \lambda using opt_tau:
  Send_Call$"tauseq" <- opt_tau
  Send_Call$"eps" <- 1e-6
  opt_lambda <- do.call(FRiSO_univar2wass, args = Send_Call)[, 1]

  # Identify selected variables at opt_tau, using given threshold
  opt_selected <- which(opt_lambda > thresh)

  # Return outputs:
  return(list(
    "tauseq" = full_tauseq,
    "errors" = errors,
    "error_sum" = error_sum,
    "opt_tau" = opt_tau,
    "opt_lambda" = opt_lambda,
    "opt_selected" = opt_selected
  ))
}

#' Fréchet Ridge Selection Operator (FRiSO)
#' 
#' @description
#' This function calculates the Fréchet Ridge Selection Operator
#' (FRiSO; \insertCite{tucker_variable_2023}{fastfrechet}), for the space of
#' univariate distribution responses equipped with the 2-Wasserstein metric. 
#' This implements with the geodesic gradient descent method of 
#' \insertCite{coulter_fast_2024}{fastfrechet}. Observed distributions are 
#' assumed to be quantile functions evaluated on a common, equally spaced 
#' `m`-grid in \eqn{(0, 1)}. Options include user-specified total penalty `tauseq`,
#' which can take a vector input (recommended decreasing sequence) where FRiSO
#' is solved at each entry with warm starts; initial sparsity vector; box
#' constraints; and gradient descent tuning parameters, including adjustable
#' impulse parameter to utilize gradient descent with momentum.
#' 
#'
#' @inheritParams frechetreg_univar2wass
#' @param tauseq A numeric vector containing 'tau' values at
#'  which to solve FRiSO problem.
#' @param lambda_init An optional `p`-long numeric vector giving the 
#'  initial allowance vector' lambda for FRiSO algorithm (default `NULL`);
#'  will be scaled to sit on tau-simplex for first entry of `tauseq`.
#' @param eps A non-negative error tolerance parameter (default `1e-5`).
#' @param nudge A non-negative numeric scalar to offset warm starts to avoid
#'  spurious boundary values (default `0`).
#' @param alpha A non-negative dampening parameter (default `0.9`).
#' @param max_iter An integer giving the maximum number of iterations for the
#'  algorithm to run (default `1000`).
#' @param max_theta A step-size scalar parameter no larger than `pi / 4`
#'  (default `pi / 4`).
#' @param impulse A scalar between `0` and `1` which controls the "impulse" in
#'  gradient descent with momentum (default `1`). `impulse` equal to `1` means
#'  no momentum and `impulse <1` means with momentum.
#'
#'
#' @references 
#' \insertRef{coulter_fast_2024}{fastfrechet}
#' 
#' \insertRef{tucker_variable_2023}{fastfrechet}
#' 
#' @return A (`p` \eqn{\times} `length(tauseq)`) matrix column-wise containing
#'  fitted 'allowance vectors' lambda per 'tau' in `tauseq`.
#' @export
#'
#' @examples
#' # Generate data for X and Y inputs by using the output of `generate_zinbinom_qf`
#' n = 100  # number of samples - nrow(X) and nrow(Y).
#' p = 10   # number of covariates - ncol(X).
#' m = 100  # EQF grid density - ncol(Y).
#' 
#' set.seed(31)
#' mydata = fastfrechet::generate_zinbinom_qf(n = n,
#'                                            p = p,
#'                                            m = m)
#'
#' X = mydata$X  # (n x p) matrix of covariates
#' Y = mydata$Y  # (n x m) matrix of EQFs, stored row-wise
#' 
#' # Dense grid of "allowance" totals:
#' tauseq = seq(0.2, 20, 0.2)
#' 
#' # Generate estimated "allowance vector"s \lambda for each \tau, stored
#' # column-wise in matrix `L`:
#' L = fastfrechet::FRiSO_univar2wass(X = X,
#'                                    Y = Y,
#'                                    lower = 0,
#'                                    upper = Inf,
#'                                    tauseq = tauseq,
#'                                    eps = 0.001,
#'                                    nudge = 0.001)
#'                                  
#' # Plot FRiSO "allowance vector" solution paths:
#' plot(x = c(), y = c(), xlim = c(0, max(tauseq)), ylim = c(0, max(L)),
#'      main = "FRiSO Selection Paths", xlab = expression(tau),
#'      ylab = expression(lambda))
#' for(k in 1:p) lines(tauseq, L[k, ], col = (k %in% 1:4) + 1, lwd = 2)
#' legend("topleft", lwd = 2, col = 1:2, bty = "n",
#'        legend = c("Unrelated variable", "Model variable"))
#' 
#' # Note that with this sample size, there is sufficient information to
#' # correctly identify the first 4 variables significantly affect the QFs.
FRiSO_univar2wass = function(X,
                             Y,
                             lower = -Inf,
                             upper = Inf,
                             tauseq,
                             lambda_init = NULL,
                             eps = 1e-5,
                             nudge = 0,
                             alpha = 0.9,
                             max_iter = 1000,
                             max_theta = pi / 4,
                             impulse = 1.0){
  
  # Dimension and compatibility checks:
  
  # Numeric matrix checks for X and Y, Z if provided:
  check_numeric(X, "matrix", finite = TRUE)
  check_numeric(Y, "matrix", finite = TRUE)
  # Check for row matching between X and Y:
  if (nrow(X) != nrow(Y)) stop("'X' and 'Y' must have the same number of rows.")
  
  # Numeric vector and constraint checks for tauseq:
  check_numeric(tauseq, "vector", finite = TRUE)
  if(min(tauseq) <= 0) stop("'tauseq' must have positive real entries.")
  
  # Numeric vector and dimension/constraint checks for lambda, if provided:
  if (!is.null(lambda_init)){
    
    check_numeric(lambda_init, "vector", finite = TRUE)
    if (length(lambda_init) != ncol(X)) stop("'lambda' must have the same length as 'X' has columns.")
    if (any(lambda_init < 0)) stop("'lambda' must have non-negative entries.")
    
  } else {
    
    # Initialize lambda_init if not provided:
    lambda_init = rep(1, ncol(X))
    
  }
  
  # Numeric scalar check and compatibility checks for lower/upper:
  check_numeric(lower, "scalar", finite = FALSE)
  check_numeric(upper, "scalar", finite = FALSE)
  if (lower >= upper) stop("'lower' must be strictly less than 'upper'.")
  
  # Numeric scalar and constraint checks for eps:
  check_numeric(eps, "scalar", finite = TRUE)
  if (eps <= 0) stop("'eps' must be strictly positive.")
  
  # Numeric scalar and constraint checks for nudge:
  check_numeric(nudge, "scalar", finite = TRUE)
  if (nudge < 0) stop("'nudge' must be non-negative.")
  
  # Numeric scalar and constraint checks for alpha:
  check_numeric(alpha, "scalar", finite = TRUE)
  if (alpha <= 0) stop("'alpha' must be strictly positive.")
  
  # Positive integer check for max_iter:
  check_wholenumber(max_iter)
  
  # Numeric scalar and constraint checks for max_theta:
  check_numeric(max_theta, "scalar", finite = FALSE)
  if (max_theta <= 0) stop("'max_theta' must be strictly positive.")

  # Numeric scalar and constraint checks for impulse:
  check_numeric(impulse, "scalar", finite = TRUE)
  if ((impulse <= 0) | (impulse > 1)) stop("'impulse' must be within (0, 1].")

  # Check for monotonicity of Y values:
  if (ncol(Y) > 1) if (min(Y[ , -1] - Y[ , -ncol(Y)]) < 0) stop("'Y' must be row-wise monotone.")
  
  # Check for box constraints on Y values:
  if ((max(Y) > upper) | (min(Y) < lower)) stop("'Y' must obey box constraints given by 'lower' and 'upper'.")
  
  
  # Center and scale inputs:
  Xc = scaleX_cpp(X)
  
  # If all entries are zero, then cannot do FRiSO (all covariates are equivalent):
  if(all(Xc == 0)) stop("'X' must have at least one non-trivial column.")
  
  # Run FRiSO:
  friso = FRiSO_GSD(X = Xc,
                    Y = Y,
                    gamma_init = sqrt(lambda_init),
                    tauseq = tauseq,
                    lower = lower,
                    upper = upper,
                    alpha = alpha,
                    nudge = nudge * tauseq,
                    eps = eps,
                    max_iter = as.integer(max_iter),
                    max_theta = max_theta,
                    impulse = impulse)
  
  return(friso$LAMBDA)
  
}

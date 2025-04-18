#' (Global) Fréchet regression for univariate distributions with 2-Wasserstein metric
#'
#' @description
#' Calculate global Fréchet regression mean objects \insertCite{petersen_frechet_2019}{fastfrechet},
#' for the space of univariate distribution responses equipped
#' with the 2-Wasserstein metric. Observed distributions are assumed
#' to be quantile functions evaluated on a common, equally spaced `m`-grid
#' in \eqn{(0, 1)}. Default is to calculate Fréchet regression means conditioned
#' on input covariate row-vectors, but can calculate Fréchet means on different
#' set of covariate vectors. Accommodates the regularization scheme from
#' \insertCite{tucker_variable_2023}{fastfrechet}.
#'
#' @param X A (`n` \eqn{\times} `p`) "input" covariate matrix with no missing, all finite entries. We recommend you **do not** include an intercept column (i.e. all `1`'s vector).
#' @param Y A (`n` \eqn{\times} `m`) matrix of observed quantile functions evaluated on a shared uniform `m`-grid in \eqn{(0, 1)} - e.g. `seq(0.5/m, 1 - 0.5/m, len = m)` - row-wise monotone non-decreasing. Entries must obey user-specified box constraints given by `lower` and `upper` parameters.
#' @param Z An optional (`z` \eqn{\times} `p`) "output" covariate matrix (default `NULL`) with no missing, all finite entries. If the rows of \eqn{(\begin{matrix} 1 & \mathbf{Z} \end{matrix})} do not belong to the column space of \eqn{(\begin{matrix} 1 & \mathbf{X} \end{matrix})^{\top}}, then there is no unique Fréchet regression solution.
#' @param C_init An optional (`n` \eqn{\times} `m + 1`) matrix (or `z` \eqn{\times} `m + 1` if `Z` is provided) of non-negative entries, specifying initial active set(s) for optimization. Active sets are identified by positive entries, row-wise if a matrix.
#' @param lambda An optional (`p` \eqn{\times} `1`) vector with non-negative entries whose sum is strictly positive.
#'  If `NULL` a non-regularized regression is performed (the default).
#' @param lower An optional numeric scalar (default `-Inf`) lower box constraint; must be strictly less than `upper`.
#' @param upper An optional numeric scalar (default `Inf`) upper box constraint; must be strictly greater than `lower`.
#' @param eps An optional numeric scalar (default `1e-10`) error tolerance; must be strictly positive.
#' @param Ztol An optional numeric scalar (default `1e-10`) error tolerance for the row-space condition on `Z`; must be strictly positive. See Details.
#'
#' @details
#' Fréchet regression generalizes Euclidean regression to the general metric
#' space setting. In the sample setting, we observe covariate-response pairs
#' \eqn{\{(\mathbf{x}_i, \mathbf{y}_i)\} \subset \mathbb{R}^p \times \Omega},
#' where \eqn{\mathbf{x}_i}'s are rows of matrix \eqn{\mathbf{X}}
#' column-centered and scaled such that
#' \eqn{\mathrm{diag}(\mathbf{X}^{\top}\mathbf{X}) = \pmb{1}}, and \eqn{\Omega}
#' is the space of univariate quantile functions equipped with the 2-Wasserstein
#' metric \eqn{d_W(\mathbf{q}, \mathbf{p}) = \lVert\mathbf{q} - \mathbf{p}\rVert_{L^2[0,1]}}.
#' (Generally we also only observe the distributions over a discrete grid, not
#' as full functions.) The conditional Fréchet mean is estimated by
#' \deqn{\widehat{\mathbf{q}}_{\oplus}(\mathbf{z}) = \underset{\mathbf{q}\in\Omega}{\mathrm{argmin}} \sum_{j=1}^n \left\{ n^{-1} + \mathbf{z}^{\top}(\mathbf{X}^{\top}\mathbf{X})^+\mathbf{x}_j\right\} d_W^2(\mathbf{y}_j, \mathbf{q}).}
#' This is the Fréchet regression problem, a weighted-mean and quadratic
#' programming problem which this function solves with a customized dual
#' active-set method inspired by \insertCite{arnstrom_dual_2022}{fastfrechet}.
#'
#' Options include box constraints on distribution support, an option to evaluate
#' quantile function estimates on a different set of covariate vectors than the
#' "input" vectors, and an option to include an initial estimate of the active
#' constraint sets. The function accommodates non-full column rank covariate
#' matrix `X`, as well as the regularization used in the variable selection
#' method of \insertCite{tucker_variable_2023}{fastfrechet} through the
#' parameter `lambda`. By default, the function will perform non-regularized
#' regression.
#'
#' In case `X` is not full column rank, and `lambda` is not provided, the
#' Fréchet mean is unique only if the rows of `(1 Z)` are in the row-space of
#' `(1 X)`. This is because the solution would otherwise depend on the choice of
#' generalized inverse of the covariance object `X'X`. Specifically, the rows of
#' `Z` have to be *affine combinations* of the rows of `X`, i.e. linear
#' combinations where the weights sum to 1. This is more restrictive than a
#' simple linear combination due to the inclusion of an explicit intercept
#' column under the hood, which is why we ask you do not include an intercept
#' column in either `X` or `Z`. The affine combination property of `X` preserves
#' the intercept column. When either `X` is full column rank or an appropriate
#' `lambda` is specified, the covariance object is either full rank or is
#' regularized, so it has an exact and unique inverse.
#'
#' This function checks the row-space condition on `Z` if it is included and if
#' `lambda` is not included. If `Z` does not meet the row-space condition, the
#' function terminates with an error.
#'
#' @return A list object with components:
#' \tabular{ll}{
#'   `Qhat` \tab returns a (`n` \eqn{\times} `m`) matrix - or (`z` \eqn{\times} `m`) matrix if `Z` was provided - row-wise containing the solutions to the Fréchet regression problem, i.e. best-fitting quantile functions bounded between `lower` and `upper`. These quantile functions are treated as evaluated on a shared `m`-grid in \eqn{(0, 1)}, e.g. `seq(0.5/m, 1 - 0.5/m, len = m)`.\cr
#'   `Lagrange_Multiplier` \tab returns a (`n` \eqn{\times} (`m + 1`)) matrix - or (`z` \eqn{\times} (`m + 1`)) matrix if `Z` was provided - row-wise containing the Lagrange multipliers associated with the underlying QP problems. \cr
#' }
#'
#' @references
#' \insertRef{arnstrom_dual_2022}{fastfrechet}
#'
#' \insertRef{petersen_frechet_2019}{fastfrechet}
#'
#' \insertRef{tucker_variable_2023}{fastfrechet}
#'
#' @export
#'
#' @examples
#' # Generate data for X and Y inputs by using the output of `generate_zinbinom_qf`
#' n <- 100 # number of samples - nrow(X) and nrow(Y).
#' p <- 10 # number of covariates - ncol(X).
#' m <- 100 # EQF grid density - ncol(Y).
#' mseq <- seq(1 / (2 * m), 1 - 1 / (2 * m), length.out = m)
#'
#' set.seed(31)
#' mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
#'
#' X <- mydata$X # (n x p) matrix of covariates
#' Y <- mydata$Y # (n x m) matrix of EQFs, stored row-wise
#'
#' # Estimate conditional QFs:
#' output <- frechetreg_univar2wass(
#'   X = X,
#'   Y = Y,
#'   Z = NULL,
#'   C_init = NULL,
#'   lambda = NULL,
#'   lower = 0,
#'   upper = Inf
#' )
#'
#' # Note: to numerical precision, these QFs are non-decreasing...
#' min(apply(output$Qhat, 1, diff))
#'
#' # ...and bounded from below by the lower bound, zero:
#' min(output$Qhat)
#'
#' # Plot the conditional QFs:
#' plot(
#'   x = c(), y = c(), xlim = c(0, 1), ylim = c(0, max(output$Qhat)),
#'   main = "Fréchet Regression QFs", xlab = "p", ylab = "quantile"
#' )
#' for (i in 1:n) lines(mseq, output$Qhat[i, ], lwd = 2)
frechetreg_univar2wass <- function(X,
                                   Y,
                                   Z = NULL,
                                   C_init = NULL,
                                   lambda = NULL,
                                   lower = -Inf,
                                   upper = Inf,
                                   eps = 1e-10,
                                   Ztol = 1e-6) {
  # Compatibility and dimension checks:

  # Numeric matrix checks for X and Y, Z if provided:
  check_numeric(X, "matrix", finite = TRUE)
  check_numeric(Y, "matrix", finite = TRUE)
  if (!is.null(Z)) check_numeric(Z, "matrix", finite = TRUE)
  # Check for row matching between X and Y:
  if (nrow(X) != nrow(Y)) stop("'X' and 'Y' must have the same number of rows.")

  # Numeric vector and dimension/constraint checks for lambda, if provided:
  if (!is.null(lambda)) {
    check_numeric(lambda, "vector", finite = TRUE)
    if (ncol(X) != length(lambda)) stop("'lambda' must have the same length as 'X' has columns.")
    if (any(lambda < 0)) stop("'lambda' must have non-negative entries.")
  }

  # Numeric scalar check and compatibility checks for lower/upper:
  check_numeric(lower, "scalar", finite = FALSE)
  check_numeric(upper, "scalar", finite = FALSE)
  if (lower >= upper) stop("'lower' must be strictly less than 'upper'.")

  # Numeric scalar and constraint checks for eps and Ztol:
  check_numeric(eps, "scalar", finite = TRUE)
  if (eps <= 0) stop("'eps' must be strictly positive.")
  check_numeric(Ztol, "scalar", finite = TRUE)
  if (Ztol <= 0) stop("'Ztol' must be strictly positive.")

  # Check for column matching between X and Z, if provided:
  if (!is.null(Z)) if (ncol(X) != ncol(Z)) stop("'X' and 'Z' must have the same number of columns.")

  # Numeric check and conversion for optional C_init:
  if (!is.null(C_init)) {
    # Check it is a numeric matrix; infinite entries (specifically +Inf) are OK
    check_numeric(C_init, "matrix", FALSE)

    # Check for row-matching with X or Z, depending on whether Z is specified
    if (is.null(Z)) if (nrow(C_init) != nrow(X)) stop("'X' and 'C_init' must have same number of rows, if 'Z' is not provided.")
    if (!is.null(Z)) if (nrow(C_init) != nrow(Z)) stop("'Z' and 'C_init', if both provided, must have same number of rows.")

    # Check ncol(C_init) = ncol(Y) + 1
    if (ncol(C_init) != (ncol(Y) + 1)) stop("'C_init' must have one more column than 'Y'.")

    # Convert with sign operator
    C_init <- sign(C_init)

    # If there are negative entries, exit with error
    if (min(C_init) < 0) stop("'C_init' must contain non-negative entries.")
  } else {
    # If not provided, initialize with zero-matrix, with appropriate nrow
    C_init <- if (is.null(Z)) matrix(0, nrow(X), ncol(Y) + 1) else matrix(0, nrow(Z), ncol(Y) + 1)
  }

  # Check for monotonicity of Y values:
  if (ncol(Y) > 1) if (min(Y[, -1] - Y[, -ncol(Y)]) < 0) stop("'Y' must be row-wise monotone non-decreasing.")

  # Check for box constraints on Y values:
  if ((max(Y) > upper) | (min(Y) < lower)) stop("'Y' must obey box constraints given by 'lower' and 'upper'.")

  # Solve for Yhat (i.e. unconstrained weighted mean):
  {
    # Check for presence of Z, the "output matrix":
    if (!is.null(Z)) {
      # Center and scale X, Z matrices:
      output <- scaleXZ_cpp(X, Z, tol = 1e-10)
      Xc <- output$Xc
      Zc <- output$Zc

      # Grab dimensions:
      n <- nrow(Xc)
      p <- ncol(Xc)
      nz <- nrow(Zc)

      # If lambda is not present (i.e. non-regularized regression)...
      if (is.null(lambda)) {
        output <- qr_proj(cbind(1, Xc), Y, cbind(1, Zc), tol = Ztol)
        if (!output$Converged) stop("The rows of matrix (1 Z) are not all in the row-space of (1 X).")

        Yhat <- output$Solution
      } else {
        # Scale Xc and Zc by sqrt(n) to simplify algebra with ridge penalty:
        Xcn <- Xc / sqrt(n)
        Zcn <- Zc / sqrt(n)

        # For simplicity in comments, let X = Xcn and Z = Zcn. The following
        # evaluates
        #
        # ( 1/n J + Z( X'X + D^{-1} )^{-1}X' )Y
        #
        # which requires inverting a matrix. When p is large in comparison to
        # n, we can rewrite the above expression
        #
        # ( 1/n J + ZDX'(XDX' + I)^{-1} )Y
        #
        # which involves inverting an (n x n) matrix. When n is large, we use
        #
        # ( 1/n J + ZB(BX'XB + I)^{-1}BX' )Y
        #
        # where B = D{1/2}, which involves inverting a (p x p) matrix.
        if (p > (1.1 * n)) {
          # Calculate DX':
          DX <- t(Xcn) * lambda

          # Calculate XDX' + I:
          G <- Xcn %*% DX
          diag(G) <- diag(G) + 1

          # Calculate Yhat:
          Yhat <- rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + (Zcn %*% DX) %*% solve(G, Y)
        } else {
          # Find b = sqrt(d), implicitly giving B = D{1/2}:
          root_lambda <- sqrt(lambda)

          # Calculate BX'XB + I:
          G <- crossprod(Xcn) * tcrossprod(root_lambda)
          diag(G) <- diag(G) + 1

          # Calculate Yhat:
          Yhat <- rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + Zcn %*% (root_lambda * solve(G, root_lambda * crossprod(Xcn, Y)))
        }
      }
    } else {
      # Center and scale X matrix:
      Xc <- scaleX_cpp(X)

      # Grab dimensions:
      n <- nrow(Xc)
      p <- ncol(Xc)
      m <- ncol(Y)

      # If lambda is not present (i.e. non-regularized regression)...
      if (is.null(lambda)) {
        # Evaluate the economical QR decomposition of (1 X) in order to perform
        # fast projection, since
        #
        # Yhat = QQ'Y
        QR <- qr_econ_getQR(cbind(1, Xc))

        Yhat <- QR$Q %*% crossprod(QR$Q, Y)
      } else {
        # Scale Xc sqrt(n) to simplify algebra with ridge penalty:
        Xcn <- Xc / sqrt(n)

        # For simplicity in comments, let X = Xcn. The following evaluates
        #
        # ( 1/n J + X( X'X + D^{-1} )^{+}X' )Y
        #
        # which can be rewritten based on the size of n vs. p. If p is larger
        # than n, a helpful equivalent representation is
        #
        # ( 1/n J + I - (XDX' + I)^{-1} )Y
        #
        # which involves inverting an (n x n) matrix. If n is larger than p,
        # a helpful equivalent representation is
        #
        # ( 1/n J + X(DX'X + I)^{-1}DX' )Y
        #
        # which involves inverting a (p x p) matrix.
        if (p > n) {
          # Calculate XDX' + I:
          G <- Xcn %*% (lambda * t(Xcn))
          diag(G) <- diag(G) + 1

          # Calculate Yhat:
          Yhat <- rep(1, n) %*% crossprod(rep(1 / n, n), Y) + Y - solve(G, Y)
        } else {
          # Calculate DX'X + I:
          S <- crossprod(Xcn)
          G <- lambda * S
          diag(G) <- diag(G) + 1

          # Calculate Yhat:
          Yhat <- rep(1, n) %*% crossprod(rep(1 / n, n), Y) + Xcn %*% solve(G, lambda * crossprod(Xcn, Y))
        }
      }
    }
  }

  # Solve for Qhat (i.e. constrained weighted mean):
  {
    # Lagrange multiplier from custom active set method:
    Eta <- monotoneQP_cpp(
      Y = Yhat,
      C_init = C_init,
      lower = lower,
      upper = upper,
      eps = eps
    )

    # Calculate Qhat from stability optimality condition:
    Qhat <- Yhat + (Eta[, -ncol(Eta)] - Eta[, -1])
  }

  # Return Qhat value:
  return(list("Qhat" = Qhat, "Lagrange_Multiplier" = Eta))
}

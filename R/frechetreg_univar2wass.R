#' (Global) Fréchet Regression for Univariate Distributions with 2-Wasserstein Metric
#'
#' @description
#' `frechetreg_univar2wass()` calculates global Fréchet regression mean
#' objects \insertCite{petersen_frechet_2019}{fastfrechet},
#' for the space of univariate distribution responses equipped
#' with the 2-Wasserstein metric. Observed distributions are assumed
#' to be quantile functions evaluated on a common, equally spaced `m`-grid
#' in \eqn{(0, 1)}. Options include user-specified output covariate matrix,
#' generalized ridge sparsity vector \insertCite{tucker_variable_2023}{fastfrechet},
#' and box constraints. The Fréchet regression problem in this context is reducible
#' to a quadratic programming problem; the workhorse of this function is a custom
#' active set method based on \insertCite{arnstrom_dual_2022}{fastfrechet}.
#'
#'
#' @param X A (`n` \eqn{\times} `p`) "input" covariate matrix with no missing, all finite entries.
#' @param Y A (`n` \eqn{\times} `m`) matrix of observed quantile functions, row-wise monotone non-decreasing.
#'  Entries must obey user-specified box constraints given by `lower` and `upper` parameters.
#' @param Z An optional (`z` \eqn{\times} `p`) "output" covariate matrix (default `NULL`) with no missing, all finite entries.
#' @param C_init An optional (`n` \eqn{\times} `m + 1`) matrix (or `z` \eqn{\times} `m + 1` if `Z` is provided) of non-negative entries, specifying initial active set(s) for optimization. Active sets are identified by positive entries, row-wise if a matrix.
#' @param lambda An optional (`p` \eqn{\times} `1`) vector with non-negative entries whose sum is strictly positive.
#'  If `NULL` a non-regularized regression is performed (the default).
#' @param lower An optional numeric scalar (default `-Inf`) lower box constraint; must be strictly less than `upper`.
#' @param upper An optional numeric scalar (default `Inf`) upper box constraint; must be strictly greater than `lower`.
#' @param eps An optional numeric scalar (default `1e-10`) error tolerance; must be strictly positive.
#'
#' @return A (`n` \eqn{\times} `m`) matrix that is the unique solution to the Fréchet Regression problem for univariate distribution responses.
#'  The solution has rows that are monotone decreasing that are bounded between the `lower` and `upper` arguments.
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
#' Generate data using the output of `generate_zinbinom_qf`
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
#' # i.e. Estimate conditional QFs:
#' Q = fastfrechet::frechetreg_univar2wass(X = X,
#'                                         Y = Y,
#'                                         Z = NULL,
#'                                         lambda = NULL,
#'                                         lower = 0,
#'                                         upper = Inf)
#'# Note: to numerical precision, these QFs are non-decreasing:
#'min(apply(Q, 1, diff))
#'# ...and bounded from below by zero, our lower rbound:
#'min(Q)
#'# Plot these conditional QFs beside the EQFs:
#'plot(x = c(), y = c(), xlim = c(0, 1), ylim = c(0, 60),
#'     main = 'Fréchet Regression QFs', xlab = 'p', ylab = 'quantile')
#'for(i in 1:n) lines(mseq, Q[i, ], col = 'red', lwd = 2)
frechetreg_univar2wass <- function(X,
                                   Y,
                                   Z = NULL,
                                   C_init = NULL,
                                   lambda = NULL,
                                   lower = -Inf,
                                   upper = Inf,
                                   eps = 1e-10) {
  
  # Compatibility and dimension checks:

  # Numeric matrix checks for X and Y, Z if provided:
  check_numeric(X, "matrix", finite = TRUE)
  check_numeric(Y, "matrix", finite = TRUE)
  if (!is.null(Z)) check_numeric(Z, "matrix", finite = TRUE)
  # Check for row matching between X and Y:
  if (nrow(X) != nrow(Y)) stop("'X' and 'Y' must have the same number of rows.")

  # Numeric vector and dimension/constraint checks for lambda, if provided:
  if (!is.null(lambda)){
    
    check_numeric(lambda, "vector", finite = TRUE)
    if (ncol(X) != length(lambda)) stop("'lambda' must have the same length as 'X' has columns.")
    if (any(lambda < 0)) stop("'lambda' must have non-negative entries.")
    
  }
  
  # Numeric scalar check and compatibility checks for lower/upper:
  check_numeric(lower, "scalar", finite = FALSE)
  check_numeric(upper, "scalar", finite = FALSE)
  if (lower >= upper) stop("'lower' must be strictly less than 'upper'.")

  # Numeric scalar and constraint checks for eps:
  check_numeric(eps, "scalar", finite = TRUE)
  if (eps <= 0) stop("'eps' must be strictly positive.")

  # Check for column matching between X and Z, if provided:
  if (!is.null(Z)) if (ncol(X) != ncol(Z)) stop("'X' and 'Z' must have the same number of columns.")
  
  # Numeric check and conversion for optional C_init:
  if(!is.null(C_init)){
    
    # Check it is a numeric matrix; infinite entries (specifically +Inf) are OK
    check_numeric(C_init, "matrix", FALSE)
    
    # Check for row-matching with X or Z, depending on whether Z is specified
    if(is.null(Z)) if(nrow(C_init) != nrow(X)) stop("'X' and 'C_init' must have same number of rows, if 'Z' is not provided.")
    if(!is.null(Z)) if(nrow(C_init) != nrow(Z)) stop("'Z' and 'C_init', if both provided, must have same number of rows.")
    
    # Check ncol(C_init) = ncol(Y) + 1
    if(ncol(C_init) != (ncol(Y) + 1)) stop("'C_init' must have one more column than 'Y'.")
    
    # Convert with sign operator
    C_init = sign(C_init)
    
    # If there are negative entries, exit with error
    if(min(C_init) < 0) stop("'C_init' must contain non-negative entries.")
    
  } else {
    
    # If not provided, initialize with zero-matrix, with appropriate nrow
    C_init = if(is.null(Z)) matrix(0, nrow(X), ncol(Y) + 1) else matrix(0, nrow(Z), ncol(Y) + 1)
    
  }

  # Check for monotonicity of Y values:
  if (ncol(Y) > 1) if (min(Y[ , -1] - Y[ , -ncol(Y)]) < 0) stop("'Y' must be row-wise monotone non-decreasing.")

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
        # For simplicity in comments, let X = Xc. The following  evaluates
        #
        # ( 1/n J + Z( X'X )^{+}X' )Y
        #
        # which requires either inverting X'X (if possible) or finding the
        # right singular vectors of X. In case X'X is not invertible, the
        # SVD of X or the SVD of X'X (both through which the right singular
        # vectors can be found) should be selected by comparing n vs. p.
        #
        # We re-write the expression of interest as
        #
        # ( 1/n J + ZM )Y,
        #
        # and try first to invert X'X. If it fails, we take an SVD route.
        if (p > (0.9 * n)) {
          # In case p out-scales n, SVD on X is faster than SVD on X'X:
          M <- tryCatch(solve(crossprod(Xc), crossprod(Xc, Y)),
            error = function(e) {
              S <- svd(Xc)
              g <- which(S$d > 1e-10)
              
              # If length(g) == 0, then X is essentially all-zeros, meaning should evaluate to all zeros:
              if(length(g) == 0) return(list(matrix(0, p, 1), matrix(0, 1, m)))
              
              return(list(S$v[ , g, drop = FALSE], crossprod(S$v[ , g, drop = FALSE], crossprod(Xc, Y)) / (S$d[g]^2)))
            }
          )

          if (!is.list(M)) {
            # Calculate Yhat:
            Yhat <- rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + Zc %*% M
          } else {
            # Calculate Yhat, preferring to first multiply Z by the right
            # singular vector matrix V[ , g] under the assumption n ~ nz:
            Yhat <- rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + (Zc %*% M[[1]]) %*% M[[2]]
          }
        } else {
          # In case n out-scales p, SVD on X'X is faster than SVD on X:
          M <- tryCatch(solve(crossprod(Xc), crossprod(Xc, Y)),
            error = function(e) {
              S <- svd(crossprod(Xc))
              g <- which(S$d > 1e-10)
              
              # If length(g) == 0, then X is essentially all-zeros, meaning should evaluate to all zeros:
              if(length(g) == 0) return(matrix(0, p, m))
              
              return(list(S$v[ , g, drop = FALSE], crossprod(S$v[ , g, drop = FALSE], crossprod(Xc, Y)) / (S$d[g])))
            }
          )

          if (!is.list(M)) {
            # Calculate Yhat:
            Yhat <- rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + Zc %*% M
          } else {
            # Calculate Yhat, preferring to first multiply the matrices from
            # M's output under the assumption n ~ nz:
            Yhat <- rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + Zc %*% (M[[1]] %*% M[[2]])
          }
        }
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
        # where B = D^{1/2}, which involves inverting a (p x p) matrix.
        if (p > (1.1 * n)) {
          # Calculate DX':
          DX <- t(Xcn) * lambda

          # Calculate XDX' + I:
          G <- Xcn %*% DX
          diag(G) <- diag(G) + 1

          # Calculate Yhat:
          Yhat <- rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + (Zcn %*% DX) %*% solve(G, Y)
        } else {
          # Find b = sqrt(d), implicitly giving B = D^{1/2}:
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
        # For simplicity in comments, let X = Xc. The following  evaluates
        #
        # ( 1/n J + X( X'X )^{+}X' )Y
        #
        # i.e. the projection of Y onto ColSpace(1, X). This requires either
        # inverting X'X (if possible) or finding the right singular vectors of
        # X.
        #
        # We re-write the expression of interest as
        #
        # ( 1/n J + ZM )Y,
        #
        # and compare p vs. n. If p is at least as large as n, we immediately
        # skip to the SVD method. If it is less than n, we try to invert X'X,
        # and if that fails then we do the SVD method.
        if (p >= n) {
          # Immediately calculate SVD of Xc = USV':
          S <- svd(Xc)
          g <- which(S$d > 1e-10)
          
          # If length(g) == 0, then X is essentially all-zeros, meaning should evaluate to all zeros:
          if (length(g) == 0){
            
            M <- 0
            
          } else {
            
            # Evaluate UU'Y, ordering the multiplication based on n > m:
            M <- if (n > m) S$u[ , g, drop = FALSE] %*% crossprod(S$u[ , g, drop = FALSE], Y) else tcrossprod(S$u[ , g, drop = FALSE]) %*% Y
            
          }

          # Calculate Yhat:
          Yhat <- rep(1, n) %*% crossprod(rep(1 / n, n), Y) + M
        } else {
          # Try to solve X(X'X)^{-1}X'Y through inversion method:
          M <- tryCatch(Xc %*% solve(crossprod(Xc), crossprod(Xc, Y)),
            error = function(e) {
              # If inversion fails, calculate SVD of Xc = USV':
              S <- svd(Xc)
              g <- which(S$d > 1e-10)
              
              # If length(g) == 0, then X is essentially all-zeros, meaning should evaluate to all zeros:
              if (length(g) == 0) return(0)
              
              # Evaluate UU'Y, ordering the multiplication based on n > m:
              if (n > m) S$u[ , g, drop = FALSE] %*% crossprod(S$u[ , g, drop = FALSE], Y) else tcrossprod(S$u[ , g, drop = FALSE]) %*% Y
            }
          )

          # Calculate Yhat:
          Yhat <- rep(1, n) %*% crossprod(rep(1 / n, n), Y) + M
        }
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
    Qhat <- Yhat + (Eta[ , -ncol(Eta)] - Eta[ , -1])
  }

  # Return Qhat value:
  return(list("Qhat" = Qhat, "Lagrange_Multiplier" = Eta))
}

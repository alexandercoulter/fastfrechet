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
#' @param tauseq A (`p` \eqn{\times} `1`) vector containing 'tau' values at
#'  which to solve FRiSO problem.
#' @param lambda_init An optional (`p` \eqn{\times} `1`) vector giving the 
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
  
  # Check for matrix inputs (X and Y; Z if provided):
  check_numeric(X, "matrix", finite = TRUE)
  check_numeric(Y, "matrix", finite = TRUE)
  # if (!(is.matrix(X) & is.matrix(Y)) | !(mode(X) == "numeric" & mode(Y) == "numeric")) stop("Y and X must be numeric matrices.")
  
  # Check for correct data type entries:
  # if (!all(is.finite(X)) | !all(is.finite(Y))) stop("Y and X must have finite numeric entries.")
  
  # Check for row matching between X and Y:
  if (nrow(X) != nrow(Y)) stop("'X' and 'Y' must have the same number of rows.")
  
  # Check for data structure and positivity of tauseq:
  check_numeric(tauseq, "vector", finite = TRUE)
  # if (!is.numeric(tauseq)) stop("'tauseq' must be a numeric vector.")
  # if (!is.vector(tauseq)) stop("'tauseq' must be a numeric vector.")
  # if (!all(is.finite(tauseq))) stop("'tauseq' must have positive real entries.")
  if(min(tauseq) <= 0) stop("'tauseq' must have positive real entries.")
  
  # Check for length, data structure, and non-negativity of sparsity vector, if provided:
  if (!is.null(lambda_init)){
    
    check_numeric(lambda_init, "vector", finite = TRUE)
    # if (!is.numeric(lambda_init)) stop("lambda must be a numeric vector.")
    # if (!is.vector(lambda_init)) stop("lambda must be a numeric vector.")
    # if (!all(is.finite(lambda_init))) stop("lambda must have finite, non-negative numeric entries.")
    if (length(lambda_init) != ncol(X)) stop("'lambda' must have the same length as 'X' has columns.")
    if (any(lambda_init < 0)) stop("'lambda' must have non-negative entries.")
    
  }
  
  # Check for box constraint compatibility:
  check_numeric(lower, "scalar", finite = FALSE)
  check_numeric(upper, "scalar", finite = FALSE)
  # if (!is.numeric(lower) | !is.numeric(upper)) stop("'lower' and 'upper' must be numeric scalars.")
  # if (!is.vector(lower) | !is.vector(upper)) stop("'lower' and 'upper' must be numeric scalars.")
  # if (length(lower) != 1 | length(upper) != 1) stop("'lower' and 'upper' must be numeric scalars.")
  if (lower >= upper) stop("'lower' must be strictly less than 'upper'.")
  
  # Check error tolerance is numeric and strictly positive:
  check_numeric(eps, "scalar", finite = TRUE)
  # f (!is.numeric(eps)) stop("'eps' must be a numeric scalar.")
  # if (!is.vector(eps)) stop("'eps' must be a numeric scalar.")
  # if (length(eps) != 1) stop("'eps' must be a numeric scalar.")
  if (eps <= 0) stop("'eps' must be strictly positive.")
  
  # Check nudge is numeric and non-negative:
  check_numeric(nudge, "scalar", finite = TRUE)
  # if (!is.numeric(nudge)) stop("'nudge' must be a numeric scalar.")
  # if (!is.vector(nudge)) stop("'nudge' must be a numeric scalar.")
  # if (length(nudge) != 1) stop("'nudge' must be a numeric scalar.")
  if (nudge < 0) stop("'nudge' must be non-negative.")
  
  # Check alpha is numeric and positive:
  check_numeric(alpha, "scalar", finite = TRUE)
  # if (!is.numeric(alpha)) stop("'alpha' must be a numeric scalar.")
  # if (!is.vector(alpha)) stop("'alpha' must be a numeric scalar.")
  # if (length(alpha) != 1) stop("'alpha' must be a numeric scalar.")
  if (alpha <= 0) stop("'alpha' must be strictly positive.")
  
  # Check max_iter is an integer (or numeric equivalent) and positive.
  check_wholenumber(max_iter)
  # if(!is.vector(max_iter)) stop("'max_iter' must be a positive integer.")
  # if(length(max_iter) != 1) stop("'max_iter' must be a positive integer.")
  # if(max_iter != as.integer(max_iter)) stop("'max_iter' must be a positive integer.")
  # if(max_iter < 1) stop("'max_iter' must be a positive integer.")
  
  # Check max_theta is numeric and positive:
  check_numeric(max_theta, "scalar", finite = FALSE)
  # if (!is.numeric(max_theta)) stop("'max_theta' must be a numeric scalar.")
  # if (!is.vector(max_theta)) stop("'max_theta' must be a numeric scalar.")
  # if (length(max_theta) != 1) stop("'max_theta' must be a numeric scalar.")
  if (max_theta <= 0) stop("'max_theta' must be strictly positive.")

  # Check impulse is numeric and strictly within (0, 1]:
  check_numeric(impulse, "scalar", finite = TRUE)
  # if (!is.numeric(impulse)) stop("'impulse' must be a numeric scalar.")
  # if (!is.vector(impulse)) stop("'impulse' must be a numeric scalar.")
  # if (length(impulse) != 1) stop("'impulse' must be a numeric scalar.")
  if ((impulse <= 0) | (impulse > 1)) stop("'impulse' must be within (0, 1].")

  # Check for monotonicity of Y values:
  if (ncol(Y) > 1) if (min(Y[ , -1] - Y[ , -ncol(Y)]) < 0) stop("'Y' must be row-wise monotone.")
  
  # Check for box constraints on Y values:
  if ((max(Y) > upper) | (min(Y) < lower)) stop("'Y' must obey box constraints given by 'lower' and 'upper'.")
  
  
  # Center and scale inputs:
  Xc = scaleX_cpp(X)
  
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

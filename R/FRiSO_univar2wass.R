#' Fréchet Ridge Selection Operator (FRiSO)
#' 
#' @description
#' This function calculates the Fréchet Ridge Selection Operator (FRiSO; Tucker et al. 2023), for the space of univariate distribution responses equipped with the 2-Wasserstein metric. This implements with the geodesic gradient descent method of Coulter et al. (2024). Observed distributions are assumed to be quantile functions evaluated on a common, equally spaced `m`-grid in (0, 1). Options include user-specified total penalty "tauseq", which can take a vector input (recommended decreasing sequence) where FRiSO is solved at each entry with warm starts; initial sparsity vector; box constraints; and gradient descent tuning parameters, including adjustable impulse parameter to utilize gradient descent with momentum.
#' 
#'
#' @inheritParams frechetreg_univar2wass
#' @param tauseq A (`p` \eqn{\times} `1`) vector containing 'tau' values at which to solve FRiSO problem.
#' @param lambda_init An optional (`p` \eqn{\times} `1`) vector giving the initial 'allowance vector' lambda for FRiSO algorithm (default `NULL`); will be scaled to sit on tau-simplex for first entry of `tauseq`.
#' @param eps A non-negative error tolerance parameter (default `1e-5`).
#' @param nudge A non-negative numeric scalar to offset warm starts to avoid spurious boundary values (default `0`).
#' @param alpha A non-negative dampening parameter (default `0.9`).
#' @param max_iter An integer giving the maximum number of iterations for the algorithm to run (default `1000`).
#' @param max_theta A step-size scalar parameter no larger than `pi / 4` (default `pi / 4`).
#' @param impulse A scalar between `0` and `1` which controls the "impulse" in gradient descent with momentum (default `1`). `impulse` equal to `1` means no momentum and `impulse <1` means with momentum.
#'
#' @return A (`p` \eqn{\times} `length(tauseq)`) matrix column-wise containing fitted 'allowance vectors' lambda per 'tau' in tauseq.
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
  
  # Grab dimensions:
  n = nrow(X)
  m = ncol(Y)
  p = ncol(X)
  
  # Compatibility checks:
  if(n != nrow(Y)) stop("\'X\' and \'Y\' should have the same number of rows.")
  if(lower > upper) stop('Lower bound should be strictly less than upper bound.')
  
  if(is.null(lambda_init)) lambda_init = rep(tauseq[1] / p, p)
  if(length(lambda_init) != p) stop("\'lambda_init\' length and ncol(X) should match.")
  if(eps < 0) stop("\'eps\' should be non-negative.")
  
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

#' Generate Random Zero-Inflated Negative Binomial Quantile Functions
#'
#' @description
#' This function generates random zero-inflated negative binomial (zinbinom)
#' quantile functions, as functions of covariates \eqn{X_1}, \eqn{X_2}, \eqn{X_3}, and \eqn{X_4}.
#' Quantile function matrix \eqn{Y} and covariate matrix \eqn{X} are generated as in Experiment B in
#' \insertCite{coulter_fast_2024}{fastfrechet}, and are evaluated on a common,
#' equally spaced `m`-grid in \eqn{(0, 1)}. Options include number of samples
#' (covariate vectors, and response quantile functions), number of total 
#' covariates (at least 4; generated iid from \eqn{N(0, 1)}), grid density `m`,
#' and baseline zinbinom parameters.
#' 
#' @param n A positive integer. Determines the number of 
#' rows in the covariate matrix.
#' @param p A positive integer greater than or equal to 4. Determines the number of
#' columns in the covariate matrix.
#' @param m A positive integer.
#' @param zero_inflation A numeric value in the range \eqn{[0,1]}.
#' @param prob A numeric value in the range \eqn{(0,1]}.
#' @param size A numeric value strictly greater than 0.
#' 
#' @references 
#' \insertRef{coulter_fast_2024}{fastfrechet}
#' 
#' @return A list object with components:
#' \tabular{ll}{
#'   `X` \tab returns a `(n` \eqn{\times} `p)` covariate matrix. \cr
#'   `Y` \tab returns an `(n` \eqn{\times} `m)` quantile function matrix. \cr
#' }
#' @import stats
#' @export
#'
#' @examples
#' n = 100  # number of samples - nrow(X) and nrow(Y).
#' p = 10   # number of covariates - ncol(X).
#' m = 100  # EQF grid density - ncol(Y).
#' mseq = seq(1 / (2 * m), 1 - 1 / (2 * m), length.out = m)
#' 
#' set.seed(31)
#' mydata = fastfrechet::generate_zinbinom_qf(n = n,
#'                                            p = p,
#'                                            m = m)
#' 
#' X = mydata$X  # (n x p) matrix of covariates
#' Y = mydata$Y  # (n x m) matrix of EQFs, stored row-wise
#' 
#' # Plot the EQFs:
#' plot(x = c(), y = c(), xlim = c(0, 1), ylim = c(0, max(Y)),
#'      main = 'Zero-Inflated Negative Binomial QFs',
#'      xlab = 'p', ylab = 'quantile')
#' for(i in 1:n) lines(mseq, Y[i, ], lwd = 2)
generate_zinbinom_qf = function(n,
                                p,
                                m,
                                zero_inflation = 0.2,
                                prob = 0.5,
                                size = 10){
  
  # Dimension and compatibility checks:
  
  # Positive integer checks for dimensions:
  check_wholenumber(n)
  check_wholenumber(p)
  check_wholenumber(m)
  # Minimum value check for p:
  if(p < 4) stop("'p' must be a positive integer greater than or equal to 4.")

  # Numeric checks for zinbinom distribution parameters:
  check_numeric(zero_inflation, "scalar", finite = TRUE)
  check_numeric(prob, "scalar", finite = TRUE)
  check_numeric(size, "scalar", finite = TRUE)
  # Constraint checks for parameters:
  if((zero_inflation < 0) | (zero_inflation > 1)) stop("'zero_inflation' must be strictly in [0, 1].")
  if((prob <= 0) | (prob > 1)) stop("'prob' must be strictly in (0, 1].")
  if(size <= 0) stop("'size' must be strictly positive.")
  
  # Generate covariate matrix:
  X = matrix(rnorm(n * p), n, p)
  
  # Control variables:
  
  # Baseline log(size) mean:
  base_logsize = log(size)
  
  # Baseline logit(prob) of success mean:
  base_logitprob = log(prob / (1 - prob))
  
  # Baseline logit(zero_inflation) mean:
  base_logitzi = log(zero_inflation / (1 - zero_inflation))
  
  # log(size) adjustment due to X2, X3 covariates:
  adj_logsize = 0.2
  
  # logit(prob) adjustment from X1 covariate:
  adj_logitprob = 0.2
  
  # logit(zero_inflation) adjustment due to X4 covariate:
  adj_logitzi = 0.4
  
  # Baseline log(size) standard deviation:
  s_logsize = 0.15
  
  # Baseline logit(prob) standrad deviation:
  s_logitprob = 0.15
  
  # Baseline logit(zero_inflation) variability:
  s_logitzi = 0.1
  
  # Sample-specific size:
  n_size = exp(rnorm(n, base_logsize + adj_logsize * (X[ , 2] + X[ , 3]), s_logsize))
  
  # Sample-specific probability of success:
  n_prob = 1 / (1 + exp(-1 * rnorm(n, base_logitprob + adj_logitprob * X[ , 1], s_logitprob)))
  
  # Sample-specific zero inflation:
  n_zi = 1 / (1 + exp(-1 * rnorm(n, base_logitzi + adj_logitzi * X[ , 4], s_logitzi)))
  
  # Generate quantile functions:
  mseq = seq(1 / (2 * m), 1 - 1 / (2 * m), length.out = m)
  mseq = rep(mseq, each = n)
  mseq = pmax((mseq - rep(n_zi, m)) / (1 - rep(n_zi, m)), 0)
  Y = matrix(qnbinom(mseq, size = n_size, prob = n_prob), n, m)
  
  # Return X matrix and Y quantile function matrix:
  return(list('X' = X, 'Y' = Y))
  
}

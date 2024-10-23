#' Title
#'
#' @param n 
#' @param p 
#' @param m 
#' @param zero_inflation 
#' @param prob 
#' @param size 
#'
#' @return
#' @export
#'
#' @examples
generate_zinbinom_qf = function(n,
                                p,
                                m,
                                zero_inflation = 0.2,
                                prob = 0.5,
                                size = 10){
  
  # Dimension and compatibility checks:
  if(p < 4) stop('Must specify at least 4 covariate columns (p).')
  if(n < 1) stop('Must specify n ≥ 1.')
  if(m < 1) stop('Must specify m ≥ 1.')
  if((zero_inflation <= 0) | (zero_inflation >= 1)) stop('\'zero_inflation\' must be strictly in (0, 1).')
  if((prob <= 0) | (prob >= 1)) stop('\'prob\' must be strictly in (0, 1).')
  if(size <= 0) stop('\'size\' must be positive.')
  
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

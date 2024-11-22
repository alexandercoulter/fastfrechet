#' Complementary Pairs Stability Selection (CPSS) for FRiSO, for Univariate Distribution Responses Equipped with the 2-Wasserstein Metric
#'
#' @param X An [n by p] matrix column-wise consisting of predictor vectors.
#' @param Y An [n by m] matrix row-wise consisting of empirical quantile functions, each evaluated on a uniformly space m-grid on (0, 1).
#' @param B Positive integer number of complementary pair sub-samples to take of original data set (default 50).
#' @param thresh A finite positive scalar selection threshold, where \eqn{\hat{lambda}_k(\tau) > }`thresh` means the \eqn{k^{th}} variable is selected, and not selected otherwise.
#' @param ... Any required or optional parameters to `FRiSO_univar2wass` function.
#'
#' @return List containing: 'tau', a vector containing the values of tau which were fed into the function call; 'selected_variables' element, a [B by 2 by length(tauseq) by p] array containing 0's and 1's identifying which variables were selected by CPSS per sub-sample; 'selected_samples' element, a [B x 2 x length(tauseq) x n] array containing 0's and 1's identifying which samples were selected in course of CPSS; a 'stability_paths' element, a [length(tauseq) x p] matrix containing stability measures for each variable (column-wise) against given tauseq.
#' @export
#'
#' @examples
FRiSO_CPSS_univar2wass = function(X,
                                  Y,
                                  B = 50,
                                  thresh = 0.0001,
                                  ...){
  
  # Grab parameters:
  Call = list(...)
  full_tauseq = Call$tauseq
  
  # Extract dimensions:
  n = nrow(X)
  halfn = floor(n/2)
  p = ncol(X)
  m = ncol(Y)
  
  # Dimension checks:
  check_wholenumber(B)
  check_numeric(thresh, "scalar", finite = TRUE)
  if(thresh <= 0) stop("'thresh' must be positive.")
  
  if(n != nrow(Y)) stop("'X' and 'Y' must have same number of rows.")
  if(n <= 1) stop("'X' and 'Y' must have at least two rows (ideally at least log2(B) rows).")
  
  # Create empty selected variables (sv) and selected samples (ss) arrays:
  sv = array(NA, c(B, 2, length(full_tauseq), p))
  ss = array(NA, c(B, 2, length(full_tauseq), n))
  
  # Loop over tau values:
  for(t in 1:length(full_tauseq)){
    
    # Extract current tau:
    Call$'tauseq' = full_tauseq[t]
    
    # Loop over splits:
    for(b in 1:B){
      
      # Split the sample into complementary sets of size floor(n/2):
      s = sample(n)
      s1 = s[1:halfn]
      s2 = s[(halfn + 1):(2 * halfn)]
      
      # Steps per sample:
      # 1. Specify sub-sampled input X and Y for the FRiSO function
      # 2. Run FRiSO problem
      # 3. Save selected samples as 0/1 binary variables for i = 1,...,n
      # 4. Save selected covariates by choosing all where L > 0.0001
      
      # First sample:
      Call$'X' = X[s1, , drop = FALSE]
      Call$'Y' = Y[s1, , drop = FALSE]
      
      L1 = do.call(FRiSO_univar2wass, args = Call)
      
      ss[b, 1, t, ] = as.numeric((1:n) %in% s1)
      sv[b, 1, t, ] = as.numeric(L1[ , 1] > thresh)
      
      # Second sample:
      Call$'X' = X[s2, , drop = FALSE]
      Call$'Y' = Y[s2, , drop = FALSE]
      
      L2 = do.call(FRiSO_univar2wass, args = Call)
      
      ss[b, 2, t, ] = as.numeric((1:n) %in% s2)
      sv[b, 2, t, ] = as.numeric(L2[ , 1] > thresh)
      
    }
    
  }
  
  # Stability paths over full_tauseq:
  stability_paths = apply(sv, c(3, 4), mean)
  
  # Return outputs: tau sequence, selected variables, selected samples, and stability paths
  return(list('tau' = full_tauseq,
              'selected_variables' = sv,
              'selected_samples' = ss,
              'stability_paths' = stability_paths))
  
}
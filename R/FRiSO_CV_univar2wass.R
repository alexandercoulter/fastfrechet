#' FRiSO_CV_univar2wass
#'
#' @param X 
#' @param Y 
#' @param K 
#' @param ... other inputs to FRiSO
#'
#' @return
#' @export
#'
#' @examples
FRiSO_CV_univar2wass = function(X, Y, K = NULL, thresh = 0.0001, ...){
  
  # Grab parameters:
  Call = list(...)
  full_tauseq = Call$tauseq
  lower = Call$'lower'
  upper = Call$'upper'
  
  # Check numeric, and grab dimensions:
  check_numeric(X, "matrix", finite = TRUE)
  check_numeric(Y, "matrix", finite = TRUE)
  n = nrow(X)
  m = ncol(Y)
  p = ncol(X)
  
  # Dimension checks:
  if(n != nrow(Y)) stop('X and Y must have same number of rows.')
  
  # If K is not specified, then we set it equal to n, i.e. LOOCV:
  if(is.null(K)){
    
    K = n
    
  } else {
    
    check_wholenumber(K)
    if(K > n) stop("'K' must be a whole number no greater than nrow(X).")
    
  }
  
  # Check numeric for thresh:
  check_numeric(thresh, "scalar", finite = TRUE)
  if(thresh <= 0) stop("'thresh' must be a positive scalar.")
  
  # Subset into K folds:
  K_sub = sample(rep(1:K, ceiling(n / K))[1:n])
  
  # Create empty CV errors array:
  errors = matrix(NA, nrow = length(Call$tauseq), ncol = K)
  
  # Loop through folds:
  for(k in 1:K){
    
    # Collect training and test data:
    Xtrain = X[K_sub != k, , drop = FALSE]
    Ytrain = Y[K_sub != k, , drop = FALSE]
    
    Xtest  = X[K_sub == k, , drop = FALSE]
    Ytest  = Y[K_sub == k, , drop = FALSE]
    
    # Train, i.e. run FRiSO:
    Call$'X' = Xtrain
    Call$'Y' = Ytrain
    L = do.call(FRiSO_univar2wass, args = Call)
    
    # Calculate errors from test subset:
    C_init = matrix(0, nrow(Xtest), p + 1)
    for(t in 1:length(full_tauseq)){
      
      # Find chosen model support:
      selected = which(L[ , t] > thresh)
      
      # If no variables are selected, i.e. 'thresh' is very large compared to
      # min('tauseq') / ncol(X), then train on all-zeros matrix:
      if(length(selected) == 0){
        
        output = frechetreg_univar2wass(X = matrix(0, nrow(Xtrain), 1),
                                        Y = Ytrain,
                                        Z = matrix(0, nrow(Xtest), 1),
                                        C_init = C_init,
                                        lambda = NULL,
                                        lower = lower,
                                        upper = upper)

      } else {
        
        output = frechetreg_univar2wass(X = Xtrain[ , selected, drop = FALSE],
                                        Y = Ytrain,
                                        Z = Xtest[ , selected, drop = FALSE],
                                        C_init = C_init,
                                        lambda = NULL,
                                        lower = lower,
                                        upper = upper)
        
      }
      
      # Save Lagrange multiplier for next warm start:
      C_init = output$'Lagrange_Multiplier'
      
      # Calculate refitted test errors:
      errors[t, k] = sum((output$'Qhat' - Ytest)^2) / m
      
    }
    
  }
  
  # Return outputs: CV errors, aggregate CV error, and optimal tau:
  return(list('errors' = errors,
              'error_sum' = rowSums(errors),
              'opt_tau' = Call$tauseq[which.min(rowSums(errors))]))
  
}
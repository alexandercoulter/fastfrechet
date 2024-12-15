#' FRiSO_CV_univar2wass
#' 
#' @description
#' 
#' 
#'
#' @inheritParams frechetreg_univar2wass
#' @param X 
#' @param Y 
#' @param K Numeric whole number no bigger than `nrow(X)`.
#' @param thresh Positive scalar that determines the selection cutoff for the `lambda` vector in the variable selection problem. 
#' @param ... other inputs to FRiSO see [fastfrechet::FRiSO_univar2wass()].
#'
#' @return
#' @export
#'
#' @examples
FRiSO_CV_univar2wass = function(X,
                                Y,
                                K = NULL,
                                thresh = 0.0001,
                                ...){
  
  # Extract call parameters:
  Call = c(as.list(environment()), list(...))
  
  # Create call list for `FRiSO_univar2wass`, matching from Call if provided:
  Send_Call = formals(FRiSO_univar2wass)
  Send_Call[names(Call)[names(Call) %in% names(Send_Call)]] = Call[names(Call) %in% names(Send_Call)]
  
  # Extract necessary parameters for `frechetreg_univar2wass`:
  full_tauseq = Send_Call$'tauseq'
  lower = Send_Call$'lower'
  upper = Send_Call$'upper'
  
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
  errors = matrix(NA, nrow = length(full_tauseq), ncol = K)
  
  # Loop through folds:
  for(k in 1:K){
    
    # Collect training and test data:
    Xtrain = X[K_sub != k, , drop = FALSE]
    Ytrain = Y[K_sub != k, , drop = FALSE]
    
    Xtest  = X[K_sub == k, , drop = FALSE]
    Ytest  = Y[K_sub == k, , drop = FALSE]
    
    # Train, i.e. run FRiSO:
    Send_Call$'X' = Xtrain
    Send_Call$'Y' = Ytrain
    L = do.call(FRiSO_univar2wass, args = Send_Call)
    
    # Calculate errors from test subset:
    C_init = matrix(0, nrow(Xtest), m + 1)
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
  
  # Calculate aggregate CV error for each \tau:
  error_sum = rowSums(errors)
  
  # Calculate \tau which minimizes aggregate CV error:
  opt_tau = full_tauseq[which.min(error_sum)[1]]
  
  # Calculate high-precision 'allowance vector' \lambda using opt_tau:
  Send_Call$'tauseq' = opt_tau
  Send_Call$'eps' = 1e-6
  opt_lambda = do.call(FRiSO_univar2wass, args = Send_Call)[ , 1]
  
  # Identify selected variables at opt_tau, using given threshold
  opt_selected = which(opt_lambda > thresh)
  
  # Return outputs:
  return(list('tauseq' = full_tauseq,
              'errors' = errors,
              'error_sum' = error_sum,
              'opt_tau' = opt_tau,
              'opt_lambda' = opt_lambda,
              'opt_selected' = opt_selected))
  
}
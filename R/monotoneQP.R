#' R Wrapper for Custom QP Solver for Nearest Empirical Quantile Function in Frobenius Norm
#' 
#' @description
#' This function wraps an Rcpp implementation of a custom active set method for solving a specific quadratic programming (QP) problem, of finding the nearest vector in Frobenius norm to an input vector, with the constraints that the output vector should be monotone non-decreasing across its entries, and obey user-specified box constraints. The function permits a matrix input, where this QP problem is solved row-wise. Implementation is based on the dual-active set method of Arnström et al. 2022, taking advantage of simplifications in this setting which avoid LDL' decompositions and other costly matrix product operations.
#'
#' @param Y 
#' @param lower 
#' @param upper 
#' @param eps tolerance level
#'
#' @return
#' @export
#'
#' @examples
monotoneQP = function(Y,
                      lower = -Inf,
                      upper = Inf,
                      eps = 1e-10){
  
  # Compatibility and dimension checks:
  {
    
    # If Y is a vector, turn into a row matrix:
    if(is.vector(Y)) Y = rbind(Y)
    
    # Check that Y is a matrix:
    if(!is.matrix(Y)) stop('\'Y\' must be a matrix.')
    
    # Check Y has at least one row and one column:
    if(prod(dim(Y)) == 0) stop('\'Y\' must have at least one row and column.')
    
    # Check for box constraint compatibility:
    if(lower > upper) stop('\'lower\' must be less than or equal to \'upper\'.')
    
  }
  
  # Run custom active set method to obtain Lagrange multiplier:
  Eta = Custom_Active_Set(Y = Y,
                          L = cbind(rep(lower, nrow(Y))),
                          U = cbind(rep(upper, nrow(Y))),
                          eps = eps)
  
  # Obtain row-monotone solution:
  Q = Y + (Eta[ , -ncol(Eta)] - Eta[ , -1])
  
  # Enforce strict box constraints against any numerical violations:
  Q[Q < lower] = lower
  Q[Q > upper] = upper
  
  # Return results:
  return(list('Lagrange_Multiplier' = Eta,
              'Solution' = Q))
  
}
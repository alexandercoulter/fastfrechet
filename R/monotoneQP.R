#' R Wrapper for Custom QP Solver for Nearest Empirical Quantile Function in Frobenius Norm
#' 
#' @description
#' `monotoneQP()` wraps an [Rcpp] implementation of a custom active set method for
#'  solving a specific quadratic programming (QP) problem, of finding the nearest
#'  vector in Frobenius norm to an input vector, with the constraints that the 
#'  output vector should be monotone non-decreasing across its entries,
#'  and obey user-specified box constraints. The function permits a matrix input,
#'  where this QP problem is solved row-wise. Implementation is based on the
#'  dual-active set method of \insertCite{arnstrom_dual_2022}{fastfrechet},
#'  taking advantage of simplifications in this setting which avoid LDL' 
#'  decompositions and other costly matrix product operations.
#'  
#' @inheritParams frechetreg_univar2wass
#' @param Y A numeric vector or matrix.
#' @param eps tolerance level (default `1e-10`).
#' 
#' @references 
#' \insertRef{arnstrom_dual_2022}{fastfrechet}
#' 
#' @return A list object with components:
#' \tabular{ll}{
#'   `Lagrange_Multiplier` \tab returns a (`n` \eqn{\times} (`m + 1`)) matrix, row-wise contains the Lagrange multipliers associated with the QP problems. \cr
#'   `Solution` \tab returns a (`n` \eqn{\times} `m`) matrix, row-wise contains the solutions to the QP problems. \cr
#' }
#' 
#' @export
#'
#' @examples
monotoneQP = function(Y,
                      lower = -Inf,
                      upper = Inf,
                      eps = 1e-10){
  
  # Compatibility and dimension checks:
  
  # If Y is a vector, turn into a row matrix:
  if(is.vector(Y)) Y = rbind(Y)
  
  # Numeric matrix check for Y:
  check_numeric(Y, "matrix", finite = TRUE)
  # Check Y has at least one row and one column:
  if(prod(dim(Y)) == 0) stop("'Y' must have at least one row and column.")
  
  # Numeric scalar check and compatibility checks for lower/upper:
  check_numeric(lower, "scalar", finite = FALSE)
  check_numeric(upper, "scalar", finite = FALSE)
  if (lower >= upper) stop("'lower' must be strictly less than 'upper'.")
  
  # Numeric scalar and constraint checks for eps:
  check_numeric(eps, "scalar", finite = FALSE)
  if (eps <= 0) stop("'eps' must be strictly positive.")
  
  # Run custom active set method to obtain Lagrange multiplier 'Eta':
  Eta = Custom_Active_Set(Y = Y,
                          L = cbind(rep(lower, nrow(Y))),
                          U = cbind(rep(upper, nrow(Y))),
                          eps = eps)
  
  # Obtain row-monotone solution Q:
  Q = Y + (Eta[ , -ncol(Eta)] - Eta[ , -1])
  
  # Enforce strict box constraints against any numerical violations:
  Q[Q < lower] = lower
  Q[Q > upper] = upper
  
  # Return results:
  return(list('Lagrange_Multiplier' = Eta,
              'Solution' = Q))
  
}
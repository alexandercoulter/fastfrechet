#' R Wrapper for Custom QP Solver for Nearest Empirical Quantile Function in Frobenius Norm
#' 
#' @description
#' `monotoneQP()` wraps an [Rcpp] implementation of a custom active set method for
#'  solving a specific quadratic programming (QP) problem, of finding the nearest
#'  vector in Frobenius norm to an input vector, with the constraints that the 
#'  output vector should be monotone non-decreasing across its entries,
#'  and obey user-specified box constraints. The function permits a matrix input,
#'  where this QP problem is solved row-wise, and also permits warm starts, in
#'  the form of specifying an initial active-constraint matrix. Implementation
#'  is based on the dual-active set method of
#'  \insertCite{arnstrom_dual_2022}{fastfrechet}, taking advantage of
#'  simplifications in this setting which avoid \eqn{LDL^{\top}} decompositions
#'  and other costly matrix product operations.
#'  
#' @inheritParams frechetreg_univar2wass
#' @param Y An (`m`)-long numeric vector (or (`n` \eqn{\times} `m`) matrix) with no missing, all finite entries.
#' @param C_init An optional (`m + 1`)-long numeric vector (or (`n` \eqn{\times} `m + 1`) matrix) of non-negative entries, specifying initial active set(s) for optimization. Active sets are identified by positive entries, row-wise if a matrix.
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
#' # Set box constraints:
#' lower = 0.5
#' upper = 1.5
#' 
#' # Generate example vector:
#' m = 100
#' set.seed(31)
#' y = rnorm(m, 2 * seq(0, 1, len = m), 0.1)
#' 
#' # Calculate monotone, box-constrained projection:
#' output = monotoneQP(y, lower = lower, upper = upper)
#' 
#' # Note: to numerical precision, these QFs are non-decreasing...
#' min(diff(output$Solution[1, ]))
#' 
#' # ...bounded from below by the lower bound...
#' min(output$Solution) - lower
#' 
#' # ...and bounded from above by the upper bound...
#' max(output$Solution) - upper
#' 
#' # Plot values of the generated vector:
#' plot(y, main = "Monotone and Box Constrained", las = 1,
#'      xlab = "Vector Entry", ylab = expression("y, "*hat("q")))
#' abline(h = c(lower, upper), lty = 2, col = "gray80")
#' 
#' # Add values of monotone projection:
#' points(output$Solution[1, ], pch = 20)
#' 
#' legend("topleft", pch = c(1, 20), bty = "n",
#'        legend = c("Unconstrained", "Monotone Constrained"))
monotoneQP = function(Y,
                      C_init = NULL,
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
  
  # Numeric check and conversion for optional C_init:
  if(!is.null(C_init)){
    
    # Can take vector input; convert to row-matrix
    if(is.vector(C_init)) C_init = rbind(C_init)
    
    # Check it is a numeric matrix; infinite entries (specifically +Inf) are OK
    check_numeric(C_init, "matrix", FALSE)
    
    # Check row and column compatibility with Y
    if(nrow(C_init) != nrow(Y)) stop("'Y' and 'C_init' must have same number of rows.")
    if(ncol(C_init) != (ncol(Y) + 1)) stop("'C_init' must have one more column than 'Y'.")
    
    # Convert with sign operator
    C_init = sign(C_init)
    
    # If there are negative entries, exit with error
    if(min(C_init) < 0) stop("'C_init' must contain non-negative entries.")
    
  } else {
    
    # If not provided, initialize with zero-matrix
    C_init = matrix(0, nrow(Y), ncol(Y) + 1)
    
  }
  
  # Numeric scalar check and compatibility checks for lower/upper:
  check_numeric(lower, "scalar", finite = FALSE)
  check_numeric(upper, "scalar", finite = FALSE)
  if (lower >= upper) stop("'lower' must be strictly less than 'upper'.")
  
  # Numeric scalar and constraint checks for eps:
  check_numeric(eps, "scalar", finite = FALSE)
  if (eps <= 0) stop("'eps' must be strictly positive.")
  
  # Run custom active set method to obtain Lagrange multiplier 'Eta':
  Eta = monotoneQP_cpp(Y = Y,
                       C_init = C_init,
                       lower = lower,
                       upper = upper,
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
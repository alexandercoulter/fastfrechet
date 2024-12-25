#' Pointwise Selection Thresholds for Complementary Pairs Stability Selection (CPSS)
#' 
#' @description
#' This function calculates pointwise selection thresholds for CPSS from
#' \insertCite{shah_variable_2013}{fastfrechet}. Empirical CPSS selection proportions
#' above these thresholds have pointwise control on the number of selected variables
#' which have "low selection probability" \eqn{\theta}. This function uses the default
#' \eqn{\theta = q/p}, where `q` is the estimated model size after variable selection.
#' 
#' @param p A positive integer giving the total number of variables.
#' @param q A numeric vector giving average model size(s) after variable selection, each entry strictly within `(0, p)`.
#' @param B A positive integer giving the number of complementary pairs splits (default `50`).
#' @param E_thr A positive numeric scalar giving the desired constraint on the number of final selected variables with "low selection probability" (default `1`).
#'
#' @details
#' CPSS obtains a final selected model by identifying variables selected a high
#' proportion of times during repeated data splits. It is a meta-procedure that
#' works independent of the baseline variable selection procedure, such as LASSO
#' or, in the case of variable selection for Fréchet regression, FRiSO. Under
#' certain conditions, \insertCite{shah_variable_2013}{fastfrechet} derived
#' thresholds for this empirical selection proportion that control the expected
#' number of variables which have "low selection probability", i.e. are disfavored
#' by the baseline procedure like FRiSO. This function inverts the problem by
#' identifying the minimum threshold necessary to control the selection procedure.
#' 
#' \insertCite{shah_variable_2013}{fastfrechet} recommend defining "low selection
#' probability" as relative to model size. A variable has low selection probability
#' if it is selected less frequently than the average relative model size, that
#' is, the number of selected variables `q` divided by the number of total
#' variables `p`. This function takes `p` and `q` as inputs (`q` may be a vector),
#' as well as the number of CPSS data splits `B` and desired control level `E_thr`.
#' Note that error control is pointwise per `q`.
#' 
#' @return A list object with components:
#' \tabular{ll}{
#'   `E_thr` \tab returns the selection constraint, as provided to the function. \cr
#'   `B` \tab returns the number of complementary pairs splits, as provided to the function. \cr
#'   `relative_model_size` \tab returns the relative model size, i.e. `q / p` from the function inputs. \cr
#'   `pointwise_thresholds` \tab returns a `p`-long vector of calculated selection thresholds. See details. \cr
#' }
#' 
#' @references 
#' \insertRef{shah_variable_2013}{fastfrechet}
#' 
#' @export
#'
#' @examples
#' # Specify number of parameters and vector of post-selection model sizes:
#' p = 30
#' q = 1:10
#' 
#' # Specify number of CPSS splits, and desired error control:
#' B = 50
#' E_thr = 1
#' 
#' # Calculate thresholds:
#' output = Shah_Samworth_thresholds(p = p, q = q, B = B, E_thr = E_thr)
#' 
#' # Plot calculated thresholds against relative model sizes:
#' plot(output$relative_model_size, output$pointwise_thresholds,
#'      type = "o", pch = 20, xlab = expression(theta),
#'      ylab = expression(pi*"("*theta*")"),
#'      main = "Shah and Samworth Thresholds", ylim = c(0, 1))
Shah_Samworth_thresholds = function(p, q, B = 50, E_thr = 1){
  
  # Dimension and compatibility checks:
  check_wholenumber(p)
  check_numeric(q, "vector", finite = TRUE)
  check_wholenumber(B)
  check_numeric(E_thr, "scalar", finite = TRUE)
  if((min(q) <= 0) | (max(q) >= p)) stop("min(q) must be strictly > 0, and max(q) must be strictly < p.")
  if(E_thr <= 0) stop("E_thr must be positive.")
  
  # Model size (relative):
  theta = q / p
  
  # Create empty vector to store minimum thresholds:
  pi_min = rep(NA, length(q))
  
  # Loop over model sizes:
  for(i in 1:length(q)){
    
    # i = 1
    # Calculate minD values using provided algorithms:
    output = minD(theta[i], B, c(-1/2, -1/4))
    
    # Find which thresholds satisfy the E constraints:
    w = which((output$minD * p) <= E_thr)
    
    # If the constraint is not satisfied, choose 1; otherwise choose smallest
    # threshold which satisfies the constraint:
    pi_min[i] = if(length(w) == 0) 1 else output$thr[w[1]]
    
  }
  
  # Return output:
  return(list("E_thr" = E_thr,
              "B" = B,
              "relative_model_size" = theta,
              "pointwise_thresholds" = pi_min))
  
}

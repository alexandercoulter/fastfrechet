#' Pointwise Selection Thresholds for Complementary Pairs Stability Selection (CPSS)
#' 
#' @description
#' This function calculates pointwise selection thresholds for CPSS from \insertCite{shah_variable_2013}{fastfrechet}. Empirical CPSS selection proportions above these thresholds have pointwise control on the number of selected variables which have "low selection probability" \eqn{\theta}. This function uses the default \eqn{\theta = q/p} where `q` is the estimated model size after variable selection.
#' 
#' @param p a positive integer giving the total number of variables
#' @param q a numeric vector giving average model size(s) after variable selection, each entry strictly within `(0, p)`
#' @param B a positive integer giving the number of complementary pairs splits (default `B = 50`)
#' @param E_thr a positive scalar giving the desired constraint on the number of selected "low-probability variables" (default `E_thr = 1`)
#'
#' @return list with the following entries: (1) `E_thr`, same as input, (2) `B`, same as input, (3) `relative_model_size`, i.e. q / p from inputs, (4) `pointwise_thresholds`, a `length(q)`-long vector of calculated selection thresholds
#' @export
#'
#' @examples
Shah_Samworth_thresholds = function(p, q, B = 50, E_thr = 1){
  
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

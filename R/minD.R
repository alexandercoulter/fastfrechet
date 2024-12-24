#' Calculate min(D...) values.
#' 
#' @description
#' See Shah and Samworth (2013), eq. 8.
#'
#' @param theta false discovery rate threshold
#' @param B number of complementary pairs splits
#' @param r r-concavity parameter (negative, e.g. -1/2, -1/4)
#'
#' @return the solution to min(D...) from Shah and Samworth (2013), eq. 8.
minD = function(theta, B, r = c(-1/2, -1/4)){
  
  # Calculate min(D...) values from Shah and Samworth (2013), Appendix A4:
  vals1 = tail_prob(B, theta^2, r[1])
  vals2 = tail_prob(2 * B, theta, r[2])
  D = pmin(c(rep(1, B), vals1$Tt), vals2$Tt)
  
  # Return threshold values and min(D...) values:
  return(list("thr" = vals2$thr,
              "minD" = D))
  
}
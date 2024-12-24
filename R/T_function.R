#' Calculate \eqn{T_t(g_{a,k})}.
#'
#' @description
#' See \insertCite{shah_variable_2013}{fastfrechet}, Appendix A4.
#' 
#' @param a a value, to optimize over
#' @param k integer value between k_min, ..., B
#' @param B number of complementary pairs splits
#' @param r r-concavity parameter (negative, e.g. -1/2, -1/4)
#' @param eta false discovery rate threshold
#' @param t candidate variable selection threshold
#'
#' @return calculated \eqn{T_t(g_{a,k})} value, eq. (17) in Appendix A4 of
#' Shah and Samworth (2013)
T_function = function(a, k, B, r, eta, t){
  
  # Calculate numerator and denominator of T{t}(g{a,k}), eq. (17) in
  # Appendix A4 of Shah and Samworth (2013):
  N = (k + 1 - B * eta) * sum((a + (0:(t - 1)))^(1/r))
  D = sum((k + 1 - (0:k)) * (a + (0:k))^(1/r))
  
  # Calculate T{t}(g{a,k}):
  1 - N / D
  
}
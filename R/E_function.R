#' Calculate stability condition for \eqn{E(f_{a,k}) = \eta} from \insertCite{shah_variable_2013}{fastfrechet}, Appendix A4.
#'
#' @param loga log(a) value, to optimize over
#' @param k integer value between k_min, ..., B
#' @param B number of complementary pairs splits
#' @param r r-concavity parameter (negative, e.g. -1/2, -1/4)
#' @param eta false discovery rate threshold, e.g. \eqn{\theta = q/p} - the target mean for \eqn{E(f_{a,k})}.
#'
#' @return logit-difference between \eqn{E(f_{a,k})} and \eqn{\eta}
E_function = function(loga, k, B, r, eta){
  
  # Calculate scaling factor for logit transformation:
  m = 2 * B / k
  
  # Calculate expected value E(f{a,k}), defined in Appendix A2 of Shah and
  # Samworth (2013):
  E = sum((0:k) / B * (exp(loga) + (0:k))^(1/r)) / sum((exp(loga) + 0:k)^(1/r))
  
  # Scale for logit transformation:
  P = E * m
  
  # Stability condition in logit-scale:
  log(P) - log(1 - P) - log(eta * m) + log(1 - eta * m)
  
}
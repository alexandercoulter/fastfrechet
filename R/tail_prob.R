#' Calculate tail probabilities for r-concave functions from from \insertCite{shah_variable_2013}{fastfrechet}
#'
#' @param B number of complementary pairs splits
#' @param eta false discovery rate threshold
#' @param r r-concavity parameter (negative, e.g. -1/2, -1/4)
#'
#' @return a list containing thresholds and calculated tail probabilities
tail_prob = function(B, eta, r){
  
  # This function implements the algorithm to find the tail probabilities
  # from Shah and Samworth (2013), Appendix A4.
  
  # Calculate minimum k value:
  kmin = ceiling(B * eta * 2) + 1
  
  # Exit if kmin >= B:
  if(kmin >= B){
    
    return(list("thr" = 1:B / B,
                "Tt" = rep(1, B)))
    
  }
  
  # Create sequence of k values:
  kseq = kmin:B
  nk = length(kseq)
  
  ## ALGORITHM STEP (a)
  aMax = 20
  a_k = rep(aMax, B)
  for(ind in kseq){
    
    # Calculate a(k) which gives E(f{a,k}) = eta; note the optimization
    # occurs in log(a) scale, so we exponentiate afterward:
    a_solution = uniroot(E_function,
                         lower = -20,
                         upper = aMax,
                         k = ind,
                         B = B,
                         r = r,
                         eta = eta,
                         tol = .Machine$double.eps^0.75)$root
    a_k[ind] = a_solution
    
    # Adjust maximum value to allow more efficient optimization:
    aMax = a_solution
    
  }
  
  # Exponentiate to give a(k) scale:
  a_k = exp(a_k)
  
  ## ALGORITHM STEPS (b), (c):
  
  # Define sequence of k values to use for steps (b), (c):
  ind = kseq
  
  # Define initial vector of optimal T{t}(f*) values:
  Tt = rep(1, B)
  
  # Loop through kmin:B - the other, leading indices will be output as 1:
  for(t in kmin:B){
    
    # Calculate T on grid of "a" and "k" values givesn by "ind":
    T_evals = sapply(ind, function(i) T_function(a = a_k[i],
                                                 k = i,
                                                 B = B,
                                                 r = r,
                                                 eta = eta,
                                                 t = t))
    
    # Find which T is largest, necessary for step (b):
    w = which.max(T_evals)
    
    if(w == 1){
      
      # If w = 1, then the overall maximum occurs in the first segment, so we
      # optimize only within this region:
      if(a_k[ind[w + 1]] == a_k[ind[w]]){
        
        output2 = list("par" = a_k[ind[w + 1]],
                       "value" = T_evals[w + 1])
        
      } else {
        
        output2 = optim(par = (a_k[ind[w]] + a_k[ind[w + 1]]) / 2,
                        fn = T_function,
                        k = ind[w],
                        B = B,
                        r = r,
                        eta = eta,
                        t = t,
                        lower = a_k[ind[w + 1]],
                        upper = a_k[ind[w]],
                        method = 'Brent',
                        control = list(fnscale = -1))
        
      }
      
      Tt[t] = output2$value
      
      # We make no changes to the lower index value since it's possible the
      # subsequent solution is also in this segment.
      
    } else if(w == length(ind)){
      
      # If w = length(ind), then the overall maximum occurs in the last segment,
      # so we optimize only within this region:
      if(a_k[ind[w]] == a_k[ind[w - 1]]){
        
        output1 = list("par" = a_k[ind[w]],
                       "value" = T_evals[w])
        
      } else {
        
        output1 = optim(par = (a_k[ind[w - 1]] + a_k[ind[w]]) / 2,
                        fn = T_function,
                        k = ind[w - 1],
                        B = B,
                        r = r,
                        eta = eta,
                        t = t,
                        lower = a_k[ind[w]],
                        upper = a_k[ind[w - 1]],
                        method = 'Brent',
                        control = list(fnscale = -1))
        
      }
      
      Tt[t] = output1$value
      
      # We update the index values to include only the last 2 now, since the
      # next optimum will now occur in only the last segment.
      ind = ind[(w - 1):w]
      
    } else {
      
      # If w is a different value, the overall maximum occurs somewhere in the
      # middle, so we optimize over the two segments immediately before and
      # immediately after the grid max:
      if(a_k[ind[w]] == a_k[ind[w - 1]]){
        
        output1 = list("par" = a_k[ind[w]],
                       "value" = T_evals[w])
        
      } else {
        
        output1 = optim(par = (a_k[ind[w - 1]] + a_k[ind[w]]) / 2,
                        fn = T_function,
                        k = ind[w - 1],
                        B = B,
                        r = r,
                        eta = eta,
                        t = t,
                        lower = a_k[ind[w]],
                        upper = a_k[ind[w - 1]],
                        method = 'Brent',
                        control = list(fnscale = -1))
        
      }
      if(a_k[ind[w + 1]] == a_k[ind[w]]){
        
        output2 = list("par" = a_k[ind[w + 1]],
                       "value" = T_evals[w + 1])
        
      } else {
        
        output2 = optim(par = (a_k[ind[w]] + a_k[ind[w + 1]]) / 2,
                        fn = T_function,
                        k = ind[w],
                        B = B,
                        r = r,
                        eta = eta,
                        t = t,
                        lower = a_k[ind[w + 1]],
                        upper = a_k[ind[w]],
                        method = 'Brent',
                        control = list(fnscale = -1))
        
      }
      
      if(output1$value > output2$value){
        
        Tt[t] = output1$value
        
        # We update the lower index value to exclude anything prior to the
        # segment in which we found the current optimum:
        ind = ind[(w - 1):length(ind)]
        
      } else {
        
        Tt[t] = output2$value
        
        # We update the lower index value to exclude anything prior to the
        # segment in which we found the current optimum:
        ind = ind[w:length(ind)]
        
      }
      
    }
    
  }
  
  # Return thresholds and T{t}(f*) values:
  return(list("thr" = 1:B / B,
              "Tt" = Tt))
  
}
GlobalObj4h=function(x,Qin,lambda,lower,upper) {
  # lambda = rep(3, p)
  # Mres = t(IndivRidgeGloWassReg_FAST(x, Qin, x, lambda)[[1]])
  Mres=IndivRidgeGloWassReg(x, Qin, x, lambda, lower = lower, upper = upper)
  return(sum((Mres-Qin)^2))
}

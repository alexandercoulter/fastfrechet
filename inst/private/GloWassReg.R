#' @title Global Wasserstein Regression
#' 
#' @description  Global Frechet regression with respect to the Wasserstein distance.
#' 
#' @param xin An n by p matrix with input measurements of the predictors.
#' @param Qin An n by m matrix with values of quantile functions of which each row holds the quantile function values on an equispaced grid on [0, 1].
#' @param xout A k by p matrix with output measurements of the predictors.
#' @param lower A scalar with the lower bound of the support of the distribution. Default is \code{NULL}.
#' @param upper A scalar with the upper bound of the support of the distribution. Default is \code{NULL}.
#' @param Rsq A logical variable indicating whether R squared would be returned. Default is FALSE.
#' @param Qgrid A numerical vector of length m holding the probability grid on [0, 1] at which the input quantile functions take values. If \code{Rsquared} is TRUE, \code{Qgrid} is needed. Default is \code{seq(1,2*m,2)/2/m}.

GloWassReg = function(xin, Qin, xout, lower=NULL, upper=NULL, Rsquared=FALSE, Qgrid=NULL){
  
  if(is.vector(xin)){
    xin = as.matrix(xin)
  }
  if(is.vector(xout)){
    xout = as.matrix(xout)
  }
  if(nrow(xin)!=nrow(Qin))
    stop("xin and Qin should have the same number of rows.")
  if(ncol(xin)!=ncol(xout))
    stop("xin and xout should have the same number of columns.")
  if(Rsquared & is.null(Qgrid)){
    warning("Qgrid is missing and taking the default value.")
  }
  
  k = nrow(xout)
  n = nrow(xin)
  m = ncol(Qin)
  xbar = colMeans(xin)
  Sigma = cov(xin) * (n-1) / n
  invSigma = solve(Sigma)
  
  # if lower & upper are neither NULL
  A = cbind(diag(m), rep(0,m)) + cbind(rep(0,m), -diag(m))
  if(!is.null(upper) & !is.null(lower)){
    b0 = c(lower, rep(0,m-1), -upper)
  }else if(!is.null(upper)){
    A = A[,-1]
    b0 = c(rep(0,m-1), -upper)
  }else if(!is.null(lower)){
    A = A[,-ncol(A)]
    b0 = c(lower,rep(0,m-1))
  }else{
    A = A[,-c(1,ncol(A))]
    b0 = rep(0,m-1)
  }
  
  Qout = sapply(1:k, function(j){
    s = 1 + t(t(xin) - xbar) %*% invSigma %*% (xout[j,] - xbar)
    s = as.vector(s)
    gx = colMeans(Qin * s)
    
    res = do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
    return(sort(res$solution))
  })
  Qout = t(Qout)
  
  if (!Rsquared) {
    return(Qout)
  } else{
    Qmean = colMeans(Qin)
    if (k == n) {
      if (isTRUE(base::all.equal(matrix(xout-xin, nrow=n), matrix(numeric(n*ncol(xin)), nrow=n), tolerance=1e-10)))
        Qin.est = Qout
    } else {
      Qin.est = sapply(1:n, function(j){
        s = 1 + t(t(xin) - xbar) %*% invSigma %*% (xin[j,] - xbar)
        s = as.vector(s)
        gx = colMeans(Qin * s)
        res = do.call(quadprog::solve.QP, list(diag(m), gx, A, b0))
        return(sort(res$solution)) #return(res$solution)
      })
      Qin.est = t(Qin.est)
    }
    Rsq = ifelse(
      is.null(Qgrid),
      1 - base::sum(t(Qin - Qin.est)^2) / base::sum((t(Qin) - Qmean)^2),
      1 - pracma::trapz(x=Qgrid, y=colSums((Qin - Qin.est)^2)) /
        pracma::trapz(x=Qgrid, y=rowSums((t(Qin) - Qmean)^2))
    )
    if(Rsq < 0) Rsq = 0
    return(list(Qout=Qout, R.squared=Rsq))
  }
}

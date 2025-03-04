#' @title Modified coordinate descent algorithm for individually penalized ridge Frechet regression
#' 
#' @description This algorithm produces the optimal lambda for individually penalized ridge 
#' Frechet regression when the output is a density 
#' 
#' @details For algorithm details, see section 1 of the supplementary materials
#' 
#' @param x: n by p matrix of predictors
#' @param Qin: An n by m matrix with values of quantile functions of which each row holds the quantile function values on an equispaced grid on [0, 1]. 
#' @param tau: a scalar constraint (sum(lambda) = tau)
#' @param lambdainit: a p by 1 vector of lambdas to initialize the coordinate descent. Default is \code{NULL}.

lambdaGlobalcoordesc=function(x,Qin, tau, lambdainit=NULL, lower = NULL, upper = NULL) {
  
  # tau = 30
  # lambda_init
  p=ncol(x)
  
  # (1) get all lambdas to sum to tau
  # If not initialized, give equal weights to all
  if(is.null(lambdainit)) {
    lambdacur=tau*rep(1,p) / p 
  } else {
    lambdacur=tau*lambdainit/sum(lambdainit)
  }
  
  # (2) Run Individual Ridge Regression with current lambda
  curobj=GlobalObj4h(x,Qin,lambdacur, lower = lower,upper = upper)
  
  wu=1
  mycount=1
  eps0=1e-4
  # time1 = Sys.time()
  # (3) Search for optimal lambdas as long as max(abs(lambdacur-lambdaold)) > 1e-4 and mycount < 20
  while(wu){
    
    objold=curobj
    lambdaold=lambdacur
    
    for (j in 1:p){
      lambdacurMj=lambdacur
      # make jth lambda 0
      lambdacurMj[j]=0
      
      # If all of the others are large still
      if(sum(lambdacurMj)>eps0){ 
        # normalize lambdas
        lambdacurMj=lambdacurMj/sum(lambdacurMj)
        
        # Make vector with only jth lambda important
        lambdaj=0*lambdacur
        lambdaj[j]=1
        
        # individually penalized ridge Frechet regression to optimized over weight ttt in (0,1)
        fffun<-function(ttt){
          # Using their function for calculating the objective:
          return(GlobalObj4h(x, Qin, ((lambdaj-lambdacurMj)*ttt+lambdacurMj)*tau, lower = lower, upper = upper))
          
          # Using my function for calculating the objective:
          # return(sum((Qin - EmbeddedProblem(X = x, Y = Qin, gam = sqrt(((lambdaj-lambdacurMj)*ttt+lambdacurMj)*tau)))^2))
        }
        
        tttmin=optimize(fffun,c(0,1),tol=0.0001)
        ttt=tttmin$minimum
        
        # If ttt large, that j important, and the others less so
        lambdacur=(lambdaj*ttt+lambdacurMj*(1-ttt))*tau
        curobj=tttmin$objective
        # print(curobj)
      } # if
    } # for
    
    # If largest difference between old and new is small enough, end while loop
    if(max(abs(lambdacur-lambdaold))<1e-4) {wu=0}
    # Or if hasn't converged after 20 iterations, end while loop
    if(mycount>100) {
      wu=0
      # print("Takes more than 20 loops to converge!!!")
    }
    # print(paste("Done with iteration ", mycount, " after ", difftime(Sys.time(), time1, units = 'mins'), " minutes."))
    mycount=mycount+1
    
  } # while
  
  return(list('lambda' = lambdacur,
              'iterations' = mycount))
}

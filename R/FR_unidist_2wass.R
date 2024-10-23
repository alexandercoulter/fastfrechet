#' (Global) Fréchet Regression for Univariate Distributions with 2-Wasserstein Metric
#' 
#' @description
#' A short description...
#' 
#'
#' @param X A `(n x p)` "input" covariate matrix with no missing, all finite entries.
#' @param Y A `(n x m)` matrix of observed quantile functions, row-wise monotone non-decreasing. Entries must obey user-specified box constraints given by `lower` and `upper` parameters.
#' @param Z An optional `(z x p)` "output" covariate matrix (default `NULL`) with no missing, all finite entries.
#' @param lambda An optional `(p x 1)` vector with non-negative entries whose sum is strictly positive.
#' @param lower An optional numeric scalar (default `-Inf`) lower box constraint; must be strictly less than `upper`.
#' @param upper An optional numeric scalar (default `Inf`) upper box constraint; must be strictly greater than `lower`.
#' @param eps An optional numeric scalar (default `1e-10`) error tolerance; must be strictly positive.
#' @param skip_checks An optional logical argument (default `FALSE`). When `TRUE` compatibility and dimension checks are skipped.
#'
#'
#' @return
#' @export
#'
#' @examples
FR_unidist_2wass = function(X,
                            Y,
                            Z = NULL,
                            lambda = NULL,
                            lower = -Inf,
                            upper = Inf,
                            eps = 1e-10,
                            skip_checks = FALSE){
  
  # Compatibility and dimension checks:
  if(!skip_checks){
    
    # Check for matrix inputs (X and Y; Z if provided):
    if(!(is.matrix(X) & is.matrix(Y))) stop('Y and X should be matrices.')
    if(!is.null(Z)) if(!is.matrix(Z)) stop('Z should be a matrix.')
    
    # Check for row matching between X and Y:
    if(nrow(X) != nrow(Y)) stop('Y and X should have the same number of rows.')
    
    # Check for length and non-negativity of sparsity vector, if provided:
    if(!is.null(lambda)) if(ncol(X) != length(lambda)) stop('lambda should have the same length as X has rows.')
    if(!is.null(lambda)) if(any(lambda < 0)) stop('All lambda entries must be non-negative.')
    
    # Check for box constraint compatibility:
    if(lower >= upper) stop('Lower bound should be strictly less than upper bound.')
    
    # Check error tolerance is strictly positive:
    if(eps <= 0) stop('Error tolerance should be strictly positive.')
    
    # Check for column matching between X and Z, if provided:
    if(!is.null(Z)) if(ncol(X) != ncol(Z)) stop('X and Z must have the same number of columns.')
    
    # Check for monotonicity of Y values:
    if(ncol(Y) > 1) if(min(Y[ , -1] - Y[ , -ncol(Y)]) < 0) stop('Y should be row-wise monotone.')
    
    # Check for box constraints on Y values:
    if((max(Y) > upper) | (min(Y) < lower)) stop('Y should obey box constraints.')
    
  }
  
  # Solve for Yhat (i.e. unconstrained weighted mean):
  {
    
    # Check for presence of Z, the "output matrix":
    if(!is.null(Z)){
      
      # Center and scale X, Z matrices:
      output = scaleXZ_cpp(X, Z, tol = 1e-10);
      Xc = output$Xc
      Zc = output$Zc
      
      n = nrow(Xc)
      p = ncol(Xc)
      nz = nrow(Zc)
      
      if(is.null(lambda)){
        
        if(p > (0.9 * n)){
          
          S = svd(Xc)
          g = S$d > 1e-10
          Yhat = rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + Zc %*% S$v[ , g] %*% (crossprod(S$v[ , g], crossprod(Xc, Y)) / (S$d[g]^2))
          
        } else {
          
          S = svd(crossprod(Xc))
          g = S$d > 1e-10
          Yhat = rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + Zc %*% S$v[ , g] %*% (crossprod(S$v[ , g], crossprod(Xc, Y)) / (S$d[g]))
          
        }
        
      } else {
        
        Xcn = Xc / sqrt(n)
        Zcn = Zc / sqrt(n)
        
        if(p > (1.1 * n)){
          
          DX = t(Xcn) * lambda
          G = Xcn %*% DX
          diag(G) = diag(G) + 1
          Yhat = rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + (Zcn %*% DX) %*% solve(G, Y)
          
        } else {
          
          root_lambda = sqrt(lambda)
          G = crossprod(Xcn) * tcrossprod(root_lambda)
          diag(G) = diag(G) + 1
          Yhat = rep(1, nz) %*% crossprod(rep(1 / n, n), Y) + Zcn %*% (root_lambda * solve(G, root_lambda * crossprod(Xcn, Y)))
          
        }
        
      }
      
    } else {
      
      Xc = scaleX_cpp(X)
      n = nrow(Xc)
      p = ncol(Xc)
      
      if(is.null(lambda)){
        
        S = svd(Xc)
        g = S$d > 1e-10
        Yhat = rep(1, n) %*% crossprod(rep(1 / n, n), Y) + S$u[ , g] %*% crossprod(S$u[ , g], Y)
        
      } else {
        
        Xcn = Xc / sqrt(n)
        
        if(p > n){
          
          G = Xcn %*% (lambda * t(Xcn))
          diag(G) = diag(G) + 1
          Yhat = rep(1, n) %*% crossprod(rep(1 / n, n), Y) + Y - solve(G, Y)
          
        } else {
          
          S = crossprod(Xcn)
          G = lambda * S
          diag(G) = diag(G) + 1
          Yhat = rep(1, n) %*% crossprod(rep(1 / n, n), Y) + Xcn %*% solve(G, lambda * crossprod(Xcn, Y))
          
        }
        
      }
      
    }
    
  }
  
  # Solve for Qhat (i.e. constrained weighted mean):
  {
    
    Eta = Custom_Active_Set(Y = Yhat,
                            L = cbind(rep(lower, nrow(Yhat))),
                            U = cbind(rep(upper, nrow(Yhat))),
                            eps = eps)
    Qhat = Yhat + (Eta[ , -ncol(Eta)] - Eta[ , -1])
    
  }
  
  # Return Qhat value:
  return(Qhat)
  
}

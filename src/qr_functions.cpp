#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List qr_econ_getQR(const arma::mat& X){
  
  arma::mat Q;
  arma::mat R;
  arma::qr_econ(Q, R, X);
  // Return:
  return(Rcpp::List::create(Rcpp::Named("Q") = Q,
                            Rcpp::Named("R") = R));
  
}

// [[Rcpp::export]]
Rcpp::List qr_proj(const arma::mat& X,
                   const arma::mat& Y,
                   const arma::mat& Z,
                   const double& tol = 1e-9){
  
  // This function finds, if it exists, a matrix X which satisfies
  //
  // Z = AX
  //
  // If such a matrix exists, then Z belongs to the column-space of X', and the
  // projection operation
  // 
  // M = Z(X'X)^{-}X'Y
  //
  // is invariant to choice of generalized inverse (X'X)^{-}. In particular, we
  // attempt to find A through the operations
  //
  // X'A' = Z'
  // QRA' = Z'
  // Q'QRA' = Q'Z'
  // RA' = Q'Z'
  // A' = R^{-1}Q'Z'
  // 
  // Note that
  // P(X)A' = X(X'X)^-X'R^{-1}Q'Z'
  //        = R'Q'(QRR'Q')^-QRR^{-1}Q'Z'
  //        = R'Q'(QRR'Q')^-Z'
  //        = R'Q'QR^{-T}R^{-1}Q'Z'
  //        = R'R^{-T}R^{-1}Q'Z'
  //        = R^{-1}Q'Z'
  //        = A',
  //
  // that is, the A' which is solved by this method is within the column space
  // of X.  If it satisfies (AX - Z) = 0, then we simply calculate
  //
  // M = Z(X'X)^{-}X'Y = AP(X)Y = AY.
  //
  // The same basic result - P(X)A' = A' - holds in the case n > p, too.
  // However, we solve instead X = QR, and then
  //
  // A' = QR^{-T}Z'
  //
  // gives the solution for A' of minimal norm (see Wikipedia).
  // 
  // Despite M = AY being easy to calculate, it can still be a computationally
  // demanding operation if nrow(X) or nrow(Z) are large.  We then look to the
  // full form of M = AY using the just-stated equation
  //
  // M = AY = ZR^{-1}Q'Y
  //
  // The matrix Q is (n x p) in dimension.  If we calculate A' first, then to
  // perform M = AY we need to perform a matrix product of two matrices with
  // dimensions (nz x n)(n x m), where we would instead like to "collapse" the
  // n dimension into a smaller r < p dimension.  We would then prefer to
  // calculate
  // 
  // M = (ZR^{-1})(Q'Y)
  //
  // because the dimensions of the two resulting matrices are (nz x p)(p x m).
  // This would be a much faster operation than a (nz x n)(n x m) operation if
  // n > p.  Luckily, in obtaining A', we do perform the left-most operation
  // 
  // B = R^{-T}Z'
  // A = B'Q'
  // 
  // We then test AX - Z = 0 and perform
  //
  // M = B'(Q'Y)
  //
  // Note that fast evaluation of AX is straightforward in the n > p case, as
  // above.  Observe
  //
  // AX = B'Q'X = B'Q'QR = B'R
  // 
  // This avoids the longer "n" dimension.
  // 
  
  int p = X.n_cols;
  int n = X.n_rows;
  
  if(p >= n){
    
    // Calculate 
    arma::mat Q(p, n);
    arma::mat R(n, n);
    arma::qr_econ(Q, R, X.t());
    arma::mat At = arma::solve(trimatu(R), (Z * Q).t());
    
    // Check the condition for AX - Z = 0:
    if((abs(At.t() * X - Z)).max() > tol){
      
      return(Rcpp::List::create(Rcpp::Named("Converged") = false,
                                Rcpp::Named("Solution") = 0,
                                Rcpp::Named("Q") = Q,
                                Rcpp::Named("R") = R));
      
    } else {
      
      return(Rcpp::List::create(Rcpp::Named("Converged") = true,
                                Rcpp::Named("Solution") = At.t() * Y,
                                Rcpp::Named("Q") = Q,
                                Rcpp::Named("R") = R));
      
    }
    
  } else {
    
    arma::mat Q(n, p);
    arma::mat R(p, p);
    arma::qr_econ(Q, R, X);
    arma::mat B = arma::solve(trimatl(R.t()), Z.t());
    
    // Check the condition for AX - Z = 0:
    if((abs(B.t() * R - Z)).max() > tol){
      
      return(Rcpp::List::create(Rcpp::Named("Converged") = false,
                                Rcpp::Named("Solution") = 0,
                                Rcpp::Named("Q") = Q,
                                Rcpp::Named("R") = R));
      
    } else {
      
      return(Rcpp::List::create(Rcpp::Named("Converged") = true,
                                Rcpp::Named("Solution") = B.t() * (Q.t() * Y),
                                Rcpp::Named("Q") = Q,
                                Rcpp::Named("R") = R));
      
    }
    
  }
  
}
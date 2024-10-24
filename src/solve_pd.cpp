#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @export
 // [[Rcpp::export]]
 arma::mat solve_pd(const arma::mat& X){
   
   return(arma::inv_sympd(X));
   
 }

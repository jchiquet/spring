#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericVector coordinate_l1(arma::vec x0, arma::vec xty, arma::mat xtx, double pen, double thr, int max_iter) {
  
  vec diag = xtx.diag()   ;

  vec xtxw = xtx * x0;
  
  vec xk = x0          ; // output vector
  int j, i = 0         ; // current iterate
  double delta = 2*thr ; // change in beta
  double u, d          ; // temporary scalar
  
  while ((delta > thr) && (i < max_iter)) {
    delta = 0;
    for (j=0; j<x0.n_elem; j++) {
      // Soft thresholding operator
      u = x0(j) * diag(j) - xty(j) - xtxw(j);
      xk(j)  = fmax(1-pen/fabs(u),0) * u/diag(j) ;
      d = xk(j)-x0(j);
      delta += pow(d,2);
      xtxw  += d*xtx.col(j) ;
    }

    // preparing next iterate
    delta = sqrt(delta);
    x0 = xk;
    i++;
    
    R_CheckUserInterrupt();
  }

  return wrap(xk) ;
  
}

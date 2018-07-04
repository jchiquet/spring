// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// coordinate_l1
Rcpp::NumericVector coordinate_l1(arma::vec x0, arma::vec xty, arma::mat xtx, double pen, double thr, int max_iter);
RcppExport SEXP _spring_coordinate_l1(SEXP x0SEXP, SEXP xtySEXP, SEXP xtxSEXP, SEXP penSEXP, SEXP thrSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xty(xtySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xtx(xtxSEXP);
    Rcpp::traits::input_parameter< double >::type pen(penSEXP);
    Rcpp::traits::input_parameter< double >::type thr(thrSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(coordinate_l1(x0, xty, xtx, pen, thr, max_iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spring_coordinate_l1", (DL_FUNC) &_spring_coordinate_l1, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_spring(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
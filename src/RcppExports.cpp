// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// comte
List comte(arma::mat y, arma::mat x, arma::mat S, double tol, int maxit, Nullable<NumericVector> min_s2, double scale, double cutoff);
RcppExport SEXP _cole_comte(SEXP ySEXP, SEXP xSEXP, SEXP SSEXP, SEXP tolSEXP, SEXP maxitSEXP, SEXP min_s2SEXP, SEXP scaleSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type min_s2(min_s2SEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(comte(y, x, S, tol, maxit, min_s2, scale, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// predict_comte
arma::mat predict_comte(List comte_obj, arma::mat newx);
RcppExport SEXP _cole_predict_comte(SEXP comte_objSEXP, SEXP newxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type comte_obj(comte_objSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type newx(newxSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_comte(comte_obj, newx));
    return rcpp_result_gen;
END_RCPP
}

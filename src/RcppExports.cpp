// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// scaleX_cpp
arma::mat scaleX_cpp(const arma::mat& X, const double& tol);
RcppExport SEXP _fastfrechet_scaleX_cpp(SEXP XSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(scaleX_cpp(X, tol));
    return rcpp_result_gen;
END_RCPP
}
// scaleXZ_cpp
Rcpp::List scaleXZ_cpp(const arma::mat& X, const arma::mat& Z, const double& tol);
RcppExport SEXP _fastfrechet_scaleXZ_cpp(SEXP XSEXP, SEXP ZSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(scaleXZ_cpp(X, Z, tol));
    return rcpp_result_gen;
END_RCPP
}
// monotoneQP_cpp
arma::mat monotoneQP_cpp(const arma::mat& Y, const arma::mat& C_init, const double& lower, const double& upper, const double& eps);
RcppExport SEXP _fastfrechet_monotoneQP_cpp(SEXP YSEXP, SEXP C_initSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C_init(C_initSEXP);
    Rcpp::traits::input_parameter< const double& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const double& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(monotoneQP_cpp(Y, C_init, lower, upper, eps));
    return rcpp_result_gen;
END_RCPP
}
// XDXt
arma::mat XDXt(const arma::mat& X, const arma::colvec& d);
RcppExport SEXP _fastfrechet_XDXt(SEXP XSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(XDXt(X, d));
    return rcpp_result_gen;
END_RCPP
}
// projA0
arma::colvec projA0(const arma::colvec& v, const arma::colvec& h);
RcppExport SEXP _fastfrechet_projA0(SEXP vSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(projA0(v, h));
    return rcpp_result_gen;
END_RCPP
}
// FRiSO_GSD
Rcpp::List FRiSO_GSD(const arma::mat& X, const arma::mat& Y, const arma::colvec& gamma_init, const arma::colvec& tauseq, const arma::colvec& nudge, const double& lower, const double& upper, const double& alpha, const double& eps, const double& max_theta, const double& impulse, const int& max_iter);
RcppExport SEXP _fastfrechet_FRiSO_GSD(SEXP XSEXP, SEXP YSEXP, SEXP gamma_initSEXP, SEXP tauseqSEXP, SEXP nudgeSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP alphaSEXP, SEXP epsSEXP, SEXP max_thetaSEXP, SEXP impulseSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type gamma_init(gamma_initSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type tauseq(tauseqSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type nudge(nudgeSEXP);
    Rcpp::traits::input_parameter< const double& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const double& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_theta(max_thetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type impulse(impulseSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(FRiSO_GSD(X, Y, gamma_init, tauseq, nudge, lower, upper, alpha, eps, max_theta, impulse, max_iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastfrechet_scaleX_cpp", (DL_FUNC) &_fastfrechet_scaleX_cpp, 2},
    {"_fastfrechet_scaleXZ_cpp", (DL_FUNC) &_fastfrechet_scaleXZ_cpp, 3},
    {"_fastfrechet_monotoneQP_cpp", (DL_FUNC) &_fastfrechet_monotoneQP_cpp, 5},
    {"_fastfrechet_XDXt", (DL_FUNC) &_fastfrechet_XDXt, 2},
    {"_fastfrechet_projA0", (DL_FUNC) &_fastfrechet_projA0, 2},
    {"_fastfrechet_FRiSO_GSD", (DL_FUNC) &_fastfrechet_FRiSO_GSD, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastfrechet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

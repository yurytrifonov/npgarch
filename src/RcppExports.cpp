// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// GARCH_lnL
List GARCH_lnL(const arma::vec par, const arma::vec y, List ind, const arma::vec arma_order, const arma::vec garch_order, std::string distribution, double degree, NumericVector knots_mult, double knots_dist, std::string return_type, Nullable<List> kernel_args, double error_value);
RcppExport SEXP _npgarch_GARCH_lnL(SEXP parSEXP, SEXP ySEXP, SEXP indSEXP, SEXP arma_orderSEXP, SEXP garch_orderSEXP, SEXP distributionSEXP, SEXP degreeSEXP, SEXP knots_multSEXP, SEXP knots_distSEXP, SEXP return_typeSEXP, SEXP kernel_argsSEXP, SEXP error_valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type ind(indSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type arma_order(arma_orderSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type garch_order(garch_orderSEXP);
    Rcpp::traits::input_parameter< std::string >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< double >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type knots_mult(knots_multSEXP);
    Rcpp::traits::input_parameter< double >::type knots_dist(knots_distSEXP);
    Rcpp::traits::input_parameter< std::string >::type return_type(return_typeSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type kernel_args(kernel_argsSEXP);
    Rcpp::traits::input_parameter< double >::type error_value(error_valueSEXP);
    rcpp_result_gen = Rcpp::wrap(GARCH_lnL(par, y, ind, arma_order, garch_order, distribution, degree, knots_mult, knots_dist, return_type, kernel_args, error_value));
    return rcpp_result_gen;
END_RCPP
}
// GARCH_aggregate
double GARCH_aggregate(const arma::vec par, const arma::vec y, List ind, const arma::vec arma_order, const arma::vec garch_order, std::string distribution, double degree, NumericVector knots_mult, double knots_dist, std::string return_type, Nullable<List> kernel_args, double error_value);
RcppExport SEXP _npgarch_GARCH_aggregate(SEXP parSEXP, SEXP ySEXP, SEXP indSEXP, SEXP arma_orderSEXP, SEXP garch_orderSEXP, SEXP distributionSEXP, SEXP degreeSEXP, SEXP knots_multSEXP, SEXP knots_distSEXP, SEXP return_typeSEXP, SEXP kernel_argsSEXP, SEXP error_valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type ind(indSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type arma_order(arma_orderSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type garch_order(garch_orderSEXP);
    Rcpp::traits::input_parameter< std::string >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< double >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type knots_mult(knots_multSEXP);
    Rcpp::traits::input_parameter< double >::type knots_dist(knots_distSEXP);
    Rcpp::traits::input_parameter< std::string >::type return_type(return_typeSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type kernel_args(kernel_argsSEXP);
    Rcpp::traits::input_parameter< double >::type error_value(error_valueSEXP);
    rcpp_result_gen = Rcpp::wrap(GARCH_aggregate(par, y, ind, arma_order, garch_order, distribution, degree, knots_mult, knots_dist, return_type, kernel_args, error_value));
    return rcpp_result_gen;
END_RCPP
}
// GARCH_lnL_grad
NumericMatrix GARCH_lnL_grad(const arma::vec par, const arma::vec y, List ind, const arma::vec arma_order, const arma::vec garch_order, std::string distribution, double degree, NumericVector knots_mult, double knots_dist, std::string return_type, Nullable<List> kernel_args, double error_value);
RcppExport SEXP _npgarch_GARCH_lnL_grad(SEXP parSEXP, SEXP ySEXP, SEXP indSEXP, SEXP arma_orderSEXP, SEXP garch_orderSEXP, SEXP distributionSEXP, SEXP degreeSEXP, SEXP knots_multSEXP, SEXP knots_distSEXP, SEXP return_typeSEXP, SEXP kernel_argsSEXP, SEXP error_valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type ind(indSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type arma_order(arma_orderSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type garch_order(garch_orderSEXP);
    Rcpp::traits::input_parameter< std::string >::type distribution(distributionSEXP);
    Rcpp::traits::input_parameter< double >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type knots_mult(knots_multSEXP);
    Rcpp::traits::input_parameter< double >::type knots_dist(knots_distSEXP);
    Rcpp::traits::input_parameter< std::string >::type return_type(return_typeSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type kernel_args(kernel_argsSEXP);
    Rcpp::traits::input_parameter< double >::type error_value(error_valueSEXP);
    rcpp_result_gen = Rcpp::wrap(GARCH_lnL_grad(par, y, ind, arma_order, garch_order, distribution, degree, knots_mult, knots_dist, return_type, kernel_args, error_value));
    return rcpp_result_gen;
END_RCPP
}
// KernelToSpline_LS
double KernelToSpline_LS(NumericVector par, NumericVector x, double degree, NumericVector knots_mult, NumericVector den_kernel, double knots_dist);
RcppExport SEXP _npgarch_KernelToSpline_LS(SEXP parSEXP, SEXP xSEXP, SEXP degreeSEXP, SEXP knots_multSEXP, SEXP den_kernelSEXP, SEXP knots_distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type knots_mult(knots_multSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type den_kernel(den_kernelSEXP);
    Rcpp::traits::input_parameter< double >::type knots_dist(knots_distSEXP);
    rcpp_result_gen = Rcpp::wrap(KernelToSpline_LS(par, x, degree, knots_mult, den_kernel, knots_dist));
    return rcpp_result_gen;
END_RCPP
}
// dkernel
NumericVector dkernel(NumericVector x, double bw, bool is_leave_one_out, double max_diff);
RcppExport SEXP _npgarch_dkernel(SEXP xSEXP, SEXP bwSEXP, SEXP is_leave_one_outSEXP, SEXP max_diffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type is_leave_one_out(is_leave_one_outSEXP);
    Rcpp::traits::input_parameter< double >::type max_diff(max_diffSEXP);
    rcpp_result_gen = Rcpp::wrap(dkernel(x, bw, is_leave_one_out, max_diff));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_npgarch_GARCH_lnL", (DL_FUNC) &_npgarch_GARCH_lnL, 12},
    {"_npgarch_GARCH_aggregate", (DL_FUNC) &_npgarch_GARCH_aggregate, 12},
    {"_npgarch_GARCH_lnL_grad", (DL_FUNC) &_npgarch_GARCH_lnL_grad, 12},
    {"_npgarch_KernelToSpline_LS", (DL_FUNC) &_npgarch_KernelToSpline_LS, 6},
    {"_npgarch_dkernel", (DL_FUNC) &_npgarch_dkernel, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_npgarch(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

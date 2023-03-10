# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

GARCH_lnL <- function(par, y, ind, arma_order, garch_order, distribution = "normal", degree = 1, knots_mult = numeric(0), knots_dist = 0.02, return_type = "aggregate", kernel_args = NULL, error_value = -1000000000000000000.0) {
    .Call(`_npgarch_GARCH_lnL`, par, y, ind, arma_order, garch_order, distribution, degree, knots_mult, knots_dist, return_type, kernel_args, error_value)
}

GARCH_aggregate <- function(par, y, ind, arma_order, garch_order, distribution = "normal", degree = 1, knots_mult = numeric(0), knots_dist = 0.02, return_type = "aggregate", kernel_args = NULL, error_value = -1000000000000000000.0) {
    .Call(`_npgarch_GARCH_aggregate`, par, y, ind, arma_order, garch_order, distribution, degree, knots_mult, knots_dist, return_type, kernel_args, error_value)
}

GARCH_lnL_grad <- function(par, y, ind, arma_order, garch_order, distribution = "normal", degree = 1, knots_mult = numeric(0), knots_dist = 0.02, return_type = "aggregate", kernel_args = NULL, error_value = -1000000000000000000.0) {
    .Call(`_npgarch_GARCH_lnL_grad`, par, y, ind, arma_order, garch_order, distribution, degree, knots_mult, knots_dist, return_type, kernel_args, error_value)
}

KernelToSpline_LS <- function(par, x, degree, knots_mult, den_kernel, knots_dist = 0.02) {
    .Call(`_npgarch_KernelToSpline_LS`, par, x, degree, knots_mult, den_kernel, knots_dist)
}

dkernel <- function(x, bw = 0.1, is_leave_one_out = FALSE, max_diff = 3.8) {
    .Call(`_npgarch_dkernel`, x, bw, is_leave_one_out, max_diff)
}


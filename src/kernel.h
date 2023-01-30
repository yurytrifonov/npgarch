#ifndef npgarch_kernel_H
#define npgarch_kernel_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;

NumericVector dkernel(NumericVector x,
                      double bw,
                      bool is_leave_one_out,
                      double max_diff);

#endif

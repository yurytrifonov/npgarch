// [[Rcpp::depends(hpa)]]
#include <RcppArmadillo.h>
#include <hpa.h>
using namespace Rcpp;

// --------------------------------------
// --------------------------------------
// --------------------------------------

// [[Rcpp::export(rng = false)]]
NumericVector dkernel(NumericVector x,
                      double bw = 0.1,
                      bool is_leave_one_out = false,
                      double max_diff = 3.8)
{
  // Currently Gaussian kernel only is provided
  
  // Prepare the value of pi number
  static const double pi = 3.14159265;

  // Perform a sorting
  NumericVector x_initial = x;
  x = clone(x);
  x.sort();
  int n_x = x.size();
  IntegerVector x_sort_ind = match(x_initial, x);
  
  double diff_lower = -max_diff;
  double diff_upper = max_diff;
  double diff = 0.0;
  
  NumericVector kde_val = NumericVector(n_x);
  
  for (int i = 0; i < n_x; i++)
  {
    if (i > 0)
    {
      for (int j = i - 1; j >= 0; j--)
      {
        diff = (x[j] - x[i]) / bw;
        if (diff >= diff_lower)
        {
          double val_tmp = std::exp(-std::pow(diff, 2) / 2);
          kde_val[i] += val_tmp;
          kde_val[j] += val_tmp;
        }
        else 
        {
          break;
        }
      }
    }
  }
  
  // Add values if is not leave one out
  if (!is_leave_one_out)
  {
    kde_val = kde_val + 1.0;
  }
  
  // Make additional adjust
  kde_val = kde_val / std::sqrt(2.0 * pi);
  kde_val = kde_val / bw;
  
  if (is_leave_one_out)
  {
    kde_val = kde_val / (n_x - 1);
  } 
  else 
  {
    kde_val = kde_val / n_x;
  }
  
  for(int i = 0; i < n_x; i++)
  {
    if (kde_val[i] < 0.000000001)
    {
      kde_val[i] = 0.000000001;
    }
  }

  // Sort in a right order
  x_sort_ind = x_sort_ind - 1;
  kde_val = kde_val[x_sort_ind];
  
  return(kde_val);
}

// [[Rcpp::depends(hpa)]]
#include <RcppArmadillo.h>
#include <hpa.h>
#include "kernel.h"
using namespace Rcpp;

// --------------------------------------
// --------------------------------------
// --------------------------------------

// [[Rcpp::export(rng = false)]]
List GARCH_lnL(const arma::vec par,
               const arma::vec y,
               List ind,
               const arma::vec arma_order,
               const arma::vec garch_order,
               std::string distribution = "normal",
               double degree = 1,
               NumericVector knots_mult = NumericVector(0),
               double knots_dist = 0.02,
               std::string return_type = "aggregate",
               Nullable<List> kernel_args = R_NilValue,
               double error_value = -1000000000000000000.0)
{
  // Store parameters based on indeсes
  int ind_mu;
  double mu;
  if (distribution != "KDE")
  {
    ind_mu = ind["mu"];
    mu = par.at(ind_mu);
  }
  
  int ind_omega;
  double omega;
  if (distribution != "KDE")
  {
    ind_omega = ind["omega"];
    omega = par.at(ind_omega);
  }
  
  arma::uvec ind_tau;
  arma::vec tau;
  if (ind.containsElementNamed("tau"))
  {
    arma::uvec ind_tmp = ind["tau"];
    ind_tau = ind_tmp;
    int n_tau = ind_tmp.size();
    if (n_tau > 0)
    {
      if (distribution == "PGN")
      {
        tau = arma::vec(n_tau + 1);
        tau.at(0) = 1;
        tau.subvec(1, n_tau) = par.elem(ind_tau);
      }
      else
      {
        tau = par.elem(ind_tau);
      }
    }
  }

  arma::uvec ind_inner_knots;
  arma::vec inner_knots;
  if (ind.containsElementNamed("inner_knots"))
  {
    arma::uvec ind_tmp = ind["inner_knots"];
    ind_inner_knots = ind_tmp;
    if (ind_inner_knots.size() > 0)
    {
      inner_knots = par.elem(ind_inner_knots);
    }
  }
  
  arma::uvec ind_ar;
  arma::vec ar;
  if (ind.containsElementNamed("ar"))
  {
    arma::uvec ind_tmp = ind["ar"];
    ind_ar = ind_tmp;
    ar = par.elem(ind_ar);
  }
  
  arma::uvec ind_ma;
  arma::vec ma;
  if (ind.containsElementNamed("ma"))
  {
    arma::uvec ind_tmp = ind["ma"];
    ind_ma = ind_tmp;
    ma = par.elem(ind_ma);
  }
  
  arma::uvec ind_alpha;
  arma::vec alpha;
  if (ind.containsElementNamed("alpha"))
  {
    arma::uvec ind_tmp = ind["alpha"];
    ind_alpha = ind_tmp;
    alpha = par.elem(ind_alpha);
  }
  
  arma::uvec ind_beta;
  arma::vec beta;
  if (ind.containsElementNamed("beta"))
  {
    arma::uvec ind_tmp = ind["beta"];
    ind_beta = ind_tmp;
    beta = par.elem(ind_beta);
  }
  
  // List to aggregate the results
  List return_list;

  // Get the number of observations
  int n = y.size();
  
  // Adjust the number of observations
  // for the lagged parts
  int n_adj_ar = n - arma_order.at(0);

  // Get the number of parameters
  int n_par = par.size();
  
  // Parameters validation
  bool is_par_incorrect = any(alpha < 0) | any(beta < 0) | (omega <= 0) |
                          ((sum(alpha) + sum(beta)) >= 0.999);
  if(is_par_incorrect)
  {
    return_list = List::create(Named("val") = error_value);
    return(return_list);
  }
  
  // Calculate the sum of beta.at(0) and alpha.at(0) values
  double alpha_sum = sum(alpha);
  double beta_sum = sum(beta);
  double alpha_beta_sum = alpha_sum + beta_sum;
  
  // Calculate the linear combination of lags
  arma::vec y_lag;
  if (arma_order.at(0) > 0)
  {
    y_lag = arma::vec(n_adj_ar);
    for (int i = 0; i < arma_order.at(0); i++)
    {
      for (int j = 0; j < n_adj_ar; j++)
      {
        y_lag.at(j) = y_lag.at(j) + 
                      ar.at(i) * y.at(j + arma_order.at(0) - i - 1);
      }
    }
  }

  // Adjust the value of the dependent variable
  // for the number of lags
  arma::vec y_adj_ar = y.subvec(arma_order.at(0), n - 1);
  
  // Estimate shocks
  arma::vec eps = y_adj_ar - mu;
  if (arma_order.at(0) > 0)
  {
    eps = eps - y_lag;
  }
  
  // Estimate unconditional volatility
  double sigma_sqr_0 = omega / (1 - alpha_beta_sum);
  double sigma_0 = sqrt(sigma_sqr_0);
  
  // ДОДЕЛАТЬ ПОТОМ
  if(distribution == "KDE")
  {
    sigma_sqr_0 = var(eps);
    sigma_0 = std::sqrt(sigma_sqr_0);
    omega = sigma_sqr_0 * (1 - alpha_beta_sum);
    mu = mean(y);
  }
  
  // Estimate conditional volatilities
  arma::vec sigma_sqr(n_adj_ar);
  sigma_sqr.at(0) = sigma_sqr_0;
  arma::vec eps_sqr = arma::pow(eps, 2);
  for (int t = 1; t < n_adj_ar; t++)
  {
    sigma_sqr.at(t) = omega;
    for (int j = 0; j < garch_order.at(0); j++)
    {
      sigma_sqr.at(t) = sigma_sqr.at(t) + alpha.at(j) * eps_sqr.at(t - 1 - j);
    }
    for (int j = 0; j < garch_order.at(1); j++)
    {
      sigma_sqr.at(t) = sigma_sqr.at(t) + beta.at(j) * sigma_sqr.at(t - 1 - j);
    }
  }
  arma::vec sigma = sqrt(sigma_sqr);

  if (return_type == "sigma")
  {
    return_list = List::create(Named("val") = sigma);
    return(return_list);
  }

  if (return_type == "eps")
  {
    return_list = List::create(Named("eps") = eps, 
                               Named("eps_std") = eps / sigma);
    return(return_list);
  }
  
  // Vectors to store likelihood information
  arma::vec L_vector(n_adj_ar);
  arma::vec lnL_vector(n_adj_ar);
    
  // Classical GARCH model
  if (distribution == "normal")
  {
    lnL_vector = arma::log_normpdf(eps / sigma) - log(sigma);
  }

  // PGN-GARCH model
  if (distribution == "PGN")
  {
    // Convert some variables
    NumericVector tau_Rcpp = wrap(tau);
    NumericVector eps_Rcpp = wrap(eps);
    NumericVector sigma_Rcpp = wrap(sigma);
    
    // Moments of PGN distribution
    NumericVector gn_exp = hpa::ehpa(NumericVector(0),
                                     tau_Rcpp, NumericVector::create(degree),
                                     NumericVector(0), NumericVector(0),
                                     NumericVector::create(0.0),
                                     NumericVector::create(1.0),
                                     NumericVector::create(1.0),
                                     false, false);
    NumericVector gn_exp2 = hpa::ehpa(NumericVector(0),
                                      tau_Rcpp, NumericVector::create(degree),
                                      NumericVector(0), NumericVector(0),
                                      NumericVector::create(0.0),
                                      NumericVector::create(1.0),
                                      NumericVector::create(2.0));
    double gn_sd = sqrt(gn_exp2[0] - std::pow(gn_exp[0], 2));

    // Standardize shocks
    NumericVector eps_std = eps_Rcpp * gn_sd / sigma_Rcpp + gn_exp[0];

    // Estimate the likelihood
    lnL_vector = hpa::dhpa(eps_std,
                           tau_Rcpp, NumericVector::create(degree),
                           NumericVector(0), NumericVector(0),
                           NumericVector::create(0.0),
                           NumericVector::create(1.0),
                           false, true, false);
    lnL_vector = lnL_vector + log(gn_sd) - log(sigma);
  }

  // SPL-GARCH model
  if (distribution == "SPL")
  {
    // Standardize shocks
    NumericVector eps_std = wrap(eps / sigma);
    // Get the number of knots and tau parameters
    int n_knots = sum(knots_mult);
    int n_knots_0 = knots_mult.size();
    int n_tau = (sum(knots_mult) - 1) - degree;

    // Assign tau parameters
    NumericVector tau_Rcpp = wrap(tau);
    // Assign knots
    NumericVector knots_0 = NumericVector(n_knots_0);
    knots_0[0] = min(eps_std);
    knots_0[n_knots_0 - 1] = max(eps_std);
    for(int i = 1; i < (n_knots_0 - 1); i++)
    {
      knots_0[i] = par.at(ind_inner_knots(i - 1));
    }

    // Assign knots considering their multiplicity
    NumericVector knots = NumericVector(n_knots);
    int knots_counter = 0;
    for (int i = 0; i < n_knots; i++)
    {
      for (int j = 0; j < knots_mult[i]; j++)
      {
        knots[knots_counter] = knots_0[i];
        knots_counter++;
      }
    }

    // Validate the distance between
    // the knots and their sorting
    NumericVector knots_ecdf = NumericVector(n_knots_0);
    knots_ecdf[0] = 1 / n;
    knots_ecdf[n_knots_0 - 1] = 1;

    for (int i = 1; i < (n_knots_0 - 1); i++)
    {
      knots_ecdf[i] = mean(knots_0[i] >= eps_std);
    }

    for (int i = 0; i < (n_knots_0 - 1); i++)
    {
      if ((knots_ecdf[i + 1] - knots_ecdf[i]) < knots_dist)
      {
        return_list = List::create(Named("val") = error_value);
        return(return_list);
      }
    }

    // Deal with b-splines
    List b = hpa::bsplineGenerate(knots, degree, false);
    List b_new = hpa::bsplineComb(b, tau_Rcpp);
    NumericMatrix b_m = b_new["m"];

    // Moments of SPL distribution
    double spl_exp = hpa::ehsa(b_m, knots, 0.0, 1.0, 1.0);
    double spl_exp2 = hpa::ehsa(b_m, knots, 0.0, 1.0, 2.0);
    double spl_sd = sqrt(spl_exp2 - std::pow(spl_exp, 2));

    // Adjust shocks
    NumericVector eps_adj = eps_std * spl_sd + spl_exp;

    // Estimate the likelihood
    lnL_vector = hpa::dhsa(eps_adj, b_m, knots, 0.0, 1.0, true);
    lnL_vector = lnL_vector + log(spl_sd) - log(sigma);
    //lnL_vector[is_infinite(lnL_vector)] = -100;
  }
  
  // Kernel based GARCH model
  if (distribution == "KDE")
  {
    double bw = tau.at(0);
    if (bw <= 0)
    {
      return_list = List::create(Named("val") = error_value);
      return(return_list);
    }

    // Standardize shocks
    NumericVector eps_std = wrap(eps / sigma);

    // Estimate kernel density using
    // leave-one-out criteria
    arma::vec kernel_den = dkernel(eps_std, bw, true, 3.8);

    // Estimate the likelihood
    lnL_vector = log(kernel_den / sigma);
  }

  // // Kernel based GARCH model
  // if (distribution == "KDE3")
  // {
  //   double bw = par[part_c + 1];
  //   par_c++;
  // 
  //   if (bw <= 0)
  //   {
  //     return_list = List::create(Named("val") = error_value);
  //     return(return_list);
  //   }
  // 
  //   // Standardize shocks
  //   NumericVector eps_std = eps / sigma;
  // 
  //   double eps_std_sd = sd(eps_std);
  //   if((eps_std_sd <= 0.99) | (eps_std_sd >= 1.01))
  //   {
  //     return_list = List::create(Named("val") = error_value);
  //     return(return_list);
  //   }
  // 
  //   // Estimate kernel density using
  //   // leave-one-out criteria
  //   NumericVector kernel_den = dkernel(eps_std, bw, true, 3.8);
  // 
  //   // Estimate the likelihood
  //   lnL_vector = log(kernel_den / sigma);
  // }
  // 
  // // Kernel based GARCH model
  // if (distribution == "KDE2")
  // {
  //   // Convert the list from Nullable
  //   List kernel_args_1;
  //   if (kernel_args.isNotNull())
  //   {
  //     List kernel_args_tmp(kernel_args);
  //     kernel_args_1 = clone(kernel_args_tmp);
  //   } else {
  //     kernel_args_1 = List::create(Named("bw") = "SJ",
  //                                  Named("kernel") = "gaussian",
  //                                  Named("n") = 512);
  //   }
  //   
  //   // Load additional environments
  //   Rcpp::Environment stats_env("package:stats");
  //   Rcpp::Function density_R = stats_env["density"];
  //   Rcpp::Function approxfun_R = stats_env["approxfun"];
  //   
  //   // Get density function arguments
  //   std::string bw_arg;
  //   if(kernel_args_1.containsElementNamed("bw"))
  //   {
  //     CharacterVector bw_arg_vec = kernel_args_1["bw"];
  //     bw_arg = bw_arg_vec[0];
  //   } else {
  //     bw_arg = "SJ";
  //     }
  //   std::string kernel_arg;
  //   if(kernel_args_1.containsElementNamed("kernel"))
  //   {
  //     CharacterVector kernel_arg_vec = kernel_args_1["kernel"];
  //     kernel_arg = kernel_arg_vec[0];
  //   } else {
  //     kernel_arg = "gaussian";
  //   }
  //   int n_arg;
  //   if(kernel_args_1.containsElementNamed("n"))
  //   {
  //     NumericVector n_arg_vec = kernel_args_1["n"];
  //     n_arg = n_arg_vec[0];
  //   } else {
  //     n_arg = 512;
  //   }
  //   
  //   // Standardize shocks
  //   NumericVector eps_std = eps / sigma;
  //   double eps_std_min = min(eps_std) * 1.01;
  //   double eps_std_max = max(eps_std) * 1.01;
  // 
  //   // Perform density estimation
  //   List kernel_List = density_R(Rcpp::_["x"] = eps_std,
  //                                Rcpp::_["bw"] = bw_arg,
  //                                Rcpp::_["kernel"] = kernel_arg,
  //                                Rcpp::_["n"] = n_arg,
  //                                Rcpp::_["from"] = eps_std_min,
  //                                Rcpp::_["to"] = eps_std_max,
  //                                Rcpp::_["na.rm"] = false);
  // 
  //   // Get density estimates
  //   Rcpp::Function kernel_fn = approxfun_R(Rcpp::_["x"] = kernel_List);
  //   NumericVector kernel_den = kernel_fn(eps_std);
  //   
  //   // Estimate the likelihood
  //   lnL_vector = log(kernel_den / sigma);
  // }
  
  if (return_type == "individual")
  {
    return_list = List::create(Named("val") = lnL_vector);
    return(return_list);
  }
  
  if (return_type == "aggregate")
  {
    NumericVector lnL_value = NumericVector::create(sum(lnL_vector));
    if (any(is_na(lnL_value) | is_nan(lnL_value)))
    {
      return_list = List::create(Named("val") = error_value);
    } 
    else 
    {
      return_list = List::create(Named("val") = lnL_value[0]);
    }
    return(return_list);
  }
  
  stop("incorrect return_type argument, it should be either sigma, individual or aggregate");
  
  return(return_list);
}

// [[Rcpp::export(rng = false)]]
double GARCH_aggregate(const arma::vec par,
                       const arma::vec y,
                       List ind,
                       const arma::vec arma_order,
                       const arma::vec garch_order,
                       std::string distribution = "normal",
                       double degree = 1,
                       NumericVector knots_mult = NumericVector(0),
                       double knots_dist = 0.02,
                       std::string return_type = "aggregate",
                       Nullable<List> kernel_args = R_NilValue,
                       double error_value = -1000000000000000000.0)
{
  List val_List = GARCH_lnL(par, y, ind, arma_order, garch_order, distribution,
                            degree, knots_mult, knots_dist, 
                            "aggregate", kernel_args, error_value);
  double val = val_List["val"];
  
  return(val);
}

// // [[Rcpp::export(rng = false)]]
// NumericVector GARCH_vectorized(NumericMatrix par,
//                                NumericVector y,
//                                std::string distribution = "normal",
//                                double degree = 1,
//                                NumericVector knots_mult = NumericVector(0),
//                                double knots_dist = 0.02,
//                                Nullable<List> kernel_args = R_NilValue)
// {
//   int n_par = par.nrow();
//   NumericVector vals = NumericVector(n_par);
// 
//   for(int i = 0; i < n_par; i++)
//   {
// 
//     NumericVector par_new = par(i, _); 
//     List val_List = GARCH_lnL(par_new, y, distribution, 
//                               degree, knots_mult, knots_dist, 
//                               "aggregate", kernel_args);
//     double val = val_List["val"];
//     vals[i] = val;
//   }
// 
//   return(vals);
// }
// 
// // --------------------------------------
// // --------------------------------------
// // --------------------------------------
// 
// [[Rcpp::export(rng = false)]]
NumericMatrix GARCH_lnL_grad(const arma::vec par,
                             const arma::vec y,
                             List ind,
                             const arma::vec arma_order,
                             const arma::vec garch_order,
                             std::string distribution = "normal",
                             double degree = 1,
                             NumericVector knots_mult = NumericVector(0),
                             double knots_dist = 0.02,
                             std::string return_type = "aggregate",
                             Nullable<List> kernel_args = R_NilValue,
                             double error_value = -1000000000000000000.0)
{
  // Store parameters based on indexes
  int ind_mu;
  double mu;
  if (distribution != "KDE")
  {
    ind_mu = ind["mu"];
    mu = par.at(ind_mu);
  }
  
  int ind_omega;
  double omega;
  if (distribution != "KDE")
  {
    ind_omega = ind["omega"];
    omega = par.at(ind_omega);
  }
  
  arma::uvec ind_tau;
  arma::vec tau;
  if (ind.containsElementNamed("tau"))
  {
    arma::uvec ind_tmp = ind["tau"];
    ind_tau = ind_tmp;
    int n_tau = ind_tmp.size();
    if (n_tau > 0)
    {
      if (distribution == "PGN")
      {
        tau = arma::vec(n_tau + 1);
        tau.at(0) = 1;
        tau.subvec(1, n_tau) = par.elem(ind_tau);
      }
      else
      {
        tau = par.elem(ind_tau);
      }
    }
  }
  
  arma::uvec ind_inner_knots;
  arma::vec inner_knots;
  if (ind.containsElementNamed("inner_knots"))
  {
    arma::uvec ind_tmp = ind["inner_knots"];
    ind_inner_knots = ind_tmp;
    if (ind_inner_knots.size() > 0)
    {
      inner_knots = par.elem(ind_inner_knots);
    }
  }
  
  arma::uvec ind_ar;
  arma::vec ar;
  if (ind.containsElementNamed("ar"))
  {
    arma::uvec ind_tmp = ind["ar"];
    ind_ar = ind_tmp;
    ar = par.elem(ind_ar);
  }
  
  arma::uvec ind_ma;
  arma::vec ma;
  if (ind.containsElementNamed("ma"))
  {
    arma::uvec ind_tmp = ind["ma"];
    ind_ma = ind_tmp;
    ma = par.elem(ind_ma);
  }
  
  arma::uvec ind_alpha;
  arma::vec alpha;
  if (ind.containsElementNamed("alpha"))
  {
    arma::uvec ind_tmp = ind["alpha"];
    ind_alpha = ind_tmp;
    alpha = par.elem(ind_alpha);
  }
  
  arma::uvec ind_beta;
  arma::vec beta;
  if (ind.containsElementNamed("beta"))
  {
    arma::uvec ind_tmp = ind["beta"];
    ind_beta = ind_tmp;
    beta = par.elem(ind_beta);
  }
  
  // List to aggregate the results
  List return_list;
  
  // Get the number of observations
  int n = y.size();
  
  // Adjust the number of observations
  // for the lagged parts
  int n_adj_ar = n - arma_order.at(0);
  
  // Get the number of parameters
  int n_par = par.size();
  
  // Parameters validation
  bool is_par_incorrect = any(alpha < 0) | any(beta < 0) | (omega <= 0) |
                          ((sum(alpha) + sum(beta)) >= 0.999);
  if(is_par_incorrect)
  {
    NumericVector vec_tmp = Rcpp::NumericVector(n_par, error_value); 
    NumericMatrix mat_tmp = Rcpp::NumericMatrix(n_par, 1, vec_tmp.begin());
    return(mat_tmp);
  }
  
  // Calculate the sum of beta.at(0) and alpha.at(0) values
  double alpha_sum = sum(alpha);
  double beta_sum = sum(beta);
  double alpha_beta_sum = alpha_sum + beta_sum;
  
  // Calculate the linear combination of lags
  arma::vec y_lag;
  if (arma_order.at(0) > 0)
  {
    y_lag = arma::vec(n_adj_ar);
    for (int i = 0; i < arma_order.at(0); i++)
    {
      for (int j = 0; j < n_adj_ar; j++)
      {
        y_lag.at(j) = y_lag.at(j) + 
          ar.at(i) * y.at(j + arma_order.at(0) - i - 1);
      }
    }
  }
  
  // Adjust the value of the dependent variable
  // for the number of lags
  arma::vec y_adj_ar = y.subvec(arma_order.at(0), n - 1);
  
  // Estimate shocks
  arma::vec eps = y_adj_ar - mu;
  if (arma_order.at(0) > 0)
  {
    eps = eps - y_lag;
  }
  
  // Estimate unconditional volatility
  double sigma_sqr_0 = omega / (1 - alpha_beta_sum);
  double sigma_0 = sqrt(sigma_sqr_0);
  
  // ДОДЕЛАТЬ ПОТОМ
  if(distribution == "KDE")
  {
    sigma_sqr_0 = var(eps);
    sigma_0 = std::sqrt(sigma_sqr_0);
    omega = sigma_sqr_0 * (1 - alpha_beta_sum);
    mu = mean(y);
  }
  
  // Estimate conditional volatilities
  arma::vec sigma_sqr(n_adj_ar);
  sigma_sqr.at(0) = sigma_sqr_0;
  arma::vec eps_sqr = arma::pow(eps, 2);
  
  for (int t = 1; t < n_adj_ar; t++)
  {
    sigma_sqr.at(t) = omega;
    for (int j = 0; j < garch_order.at(0); j++)
    {
      sigma_sqr.at(t) = sigma_sqr.at(t) + alpha.at(j) * eps_sqr.at(t - 1 - j);
    }
    for (int j = 0; j < garch_order.at(1); j++)
    {
      sigma_sqr.at(t) = sigma_sqr.at(t) + beta.at(j) * sigma_sqr.at(t - 1 - j);
    }
  }
  arma::vec sigma = sqrt(sigma_sqr);
  
  // Calculate standardized shocks
  NumericVector eps_std = wrap(eps / sigma);

  // ЗАМЕНИТЬ n на n_adj_ar в общем случае!!!!!!!!!!!!!!!!!!!!
  
  // Vectors to store Jacobian information
  NumericMatrix lnL_jacobian = NumericMatrix(n, n_par);
  
  // Derivatives of sigma

    // respect to mu
  NumericVector d_sigma_mu = NumericVector(n);
  d_sigma_mu[0] = 0;
  for(int t = 1; t < n; t++)
  {
    d_sigma_mu[t] = (beta.at(0) * sigma[t - 1] * d_sigma_mu[t - 1] -
                     alpha.at(0) * eps[t - 1]) / sigma[t];
  }

    // respect to omega
  NumericVector d_sigma_omega = NumericVector(n);
  d_sigma_omega[0] = 1 / (2 * sigma_0 * (1 - alpha.at(0) - beta.at(0)));
  for(int t = 1; t < n; t++)
  {
    d_sigma_omega[t] = (1 + beta.at(0) * 2 * sigma[t - 1] * d_sigma_omega[t - 1]) /
                       (2 * sigma[t]);
  }

    // respect to alpha.at(0)
  NumericVector d_sigma_alpha = NumericVector(n);
  d_sigma_alpha[0] = omega / (2 * sigma_0 * std::pow((1 - alpha.at(0) - beta.at(0)), 2.0));
  for(int t = 1; t < n; t++)
  {
    d_sigma_alpha[t] = (eps_sqr[t - 1] +
                        beta.at(0) * 2 * sigma[t - 1] * d_sigma_alpha[t - 1]) /
                       (2 * sigma[t]);
  }

    // respect to beta.at(0)
  NumericVector d_sigma_beta = NumericVector(n);
  d_sigma_beta[0] = d_sigma_alpha[0];
  for(int t = 1; t < n; t++)
  {
    d_sigma_beta[t] = (1 / (2 * sigma[t])) *
                       (beta.at(0) * 2 * sigma[t - 1] * d_sigma_beta[t - 1] +
                        sigma_sqr[t - 1]);
  }

  // Derivatives of standardized shocks
  NumericVector sigma_Rcpp = wrap(sigma);
    // respect to mu
  NumericVector d_eps_std_mu = -(1.0 + eps_std * d_sigma_mu) / sigma_Rcpp;

    // respect to omega
    NumericVector d_eps_std_omega = -(eps_std / sigma_Rcpp) * d_sigma_omega;

    // respect to alpha.at(0)
    NumericVector d_eps_std_alpha = -(eps_std / sigma_Rcpp) * d_sigma_alpha;

    // respect to beta.at(0)
    NumericVector d_eps_std_beta = -(eps_std / sigma_Rcpp) * d_sigma_beta;

  // Classical GARCH model
  if(distribution == "normal")
  {
    // Estimate the Jacobian
    lnL_jacobian(_, 0) = (-eps_std) * d_eps_std_mu - d_sigma_mu / sigma_Rcpp;
    lnL_jacobian(_, 1) = (-eps_std) * d_eps_std_omega - d_sigma_omega / sigma_Rcpp;
    lnL_jacobian(_, 2) = (-eps_std) * d_eps_std_alpha - d_sigma_alpha / sigma_Rcpp;
    lnL_jacobian(_, 3) = (-eps_std) * d_eps_std_beta - d_sigma_beta / sigma_Rcpp;
  }

  // PGN-GARCH model
  if(distribution == "PGN")
  {
    NumericVector tau = NumericVector(degree + 1);
    tau[0] = 1;
    for(int i = 1; i <= degree; i++)
    {
      tau[i] = par[3 + i];
    }

    // Moments of PGN distribution
    NumericVector gn_exp = hpa::ehpa(NumericVector(0),
                                     tau, NumericVector::create(degree),
                                     NumericVector(0), NumericVector(0),
                                     NumericVector::create(0.0),
                                     NumericVector::create(1.0),
                                     NumericVector::create(1.0),
                                     false, false);

    NumericVector gn_exp2 = hpa::ehpa(NumericVector(0),
                                      tau, NumericVector::create(degree),
                                      NumericVector(0), NumericVector(0),
                                      NumericVector::create(0.0),
                                      NumericVector::create(1.0),
                                      NumericVector::create(2.0));

    double gn_sd = sqrt(gn_exp2[0] - std::pow(gn_exp[0], 2));

    // Differentiate moments respect to tau (excluding tau0)
    NumericVector d_gn_exp = hpa::ehpaDiff(NumericVector(0),
                                           tau, NumericVector::create(degree),
                                           NumericVector(0), NumericVector(0),
                                           NumericVector::create(0.0),
                                           NumericVector::create(1.0),
                                           NumericVector::create(1.0),
                                           "pol_coefficients",
                                           false, false, false);
    d_gn_exp.erase(0);

    NumericVector d_gn_exp2 = hpa::ehpaDiff(NumericVector(0),
                                            tau, NumericVector::create(degree),
                                            NumericVector(0), NumericVector(0),
                                            NumericVector::create(0.0),
                                            NumericVector::create(1.0),
                                            NumericVector::create(2.0),
                                            "pol_coefficients",
                                            false, false, false);
    d_gn_exp2.erase(0);

    NumericVector d_gn_sd = (d_gn_exp2 - 2 * gn_exp[0] * d_gn_exp) / (2 * gn_sd);

    // Adjust standardize shocks
    NumericVector eps_std_adj = eps_std * gn_sd + gn_exp[0];

    // Differentiate adjusted standardized shocks
      // respect to tau
    NumericMatrix d_eps_std_adj_tau = NumericMatrix(n, degree);
    for (int i = 0; i < degree; i++)
    {
      d_eps_std_adj_tau(_, i) = (eps_std * d_gn_sd[i]) + d_gn_exp[i];
    }

    // Differentiate the log-density
      // respect to the argument
    NumericVector d_log_density_arg = hpa::dhpaDiff(eps_std_adj,
                                                    tau, NumericVector::create(degree),
                                                    NumericVector(0), NumericVector(0),
                                                    NumericVector::create(0.0),
                                                    NumericVector::create(1.0),
                                                    "x",
                                                    false, true, false);
      // respect to the tau
    NumericMatrix d_log_density_tau = hpa::dhpaDiff(eps_std_adj,
                                                    tau, NumericVector::create(degree),
                                                    NumericVector(0), NumericVector(0),
                                                    NumericVector::create(0.0),
                                                    NumericVector::create(1.0),
                                                    "pol_coefficients",
                                                    false, true, false);

    // Estimate the Jacobian
    lnL_jacobian(_, 0) = d_log_density_arg * (d_eps_std_mu * gn_sd) -
                         d_sigma_mu / sigma_Rcpp;
    lnL_jacobian(_, 1) = d_log_density_arg * (d_eps_std_omega * gn_sd) -
                         d_sigma_omega / sigma_Rcpp;
    lnL_jacobian(_, 2) = d_log_density_arg * (d_eps_std_alpha * gn_sd) -
                         d_sigma_alpha / sigma_Rcpp;
    lnL_jacobian(_, 3) = d_log_density_arg * (d_eps_std_beta * gn_sd) -
                         d_sigma_beta / sigma_Rcpp;

    for (int i = 0; i < degree; i++)
    {
      lnL_jacobian(_, i + 4) = d_log_density_tau(_, i + 1) +
                               d_log_density_arg * d_eps_std_adj_tau(_, i) +
                               d_gn_sd[i] / gn_sd;

    }

  }

  if (return_type == "individual")
  {
    return(lnL_jacobian);
  }

  NumericVector lnL_grad_vec = colSums(lnL_jacobian);
  NumericMatrix lnL_grad = NumericMatrix(n_par, 1);
  lnL_grad(_, 0) = lnL_grad_vec;

  return(lnL_grad);
}
// 
// // --------------------------------------
// // --------------------------------------
// // --------------------------------------
  
// [[Rcpp::export(rng = false)]]
double KernelToSpline_LS(NumericVector par,
                         NumericVector x,
                         double degree,
                         NumericVector knots_mult,
                         NumericVector den_kernel,
                         double knots_dist = 0.02)
{
  // Get the number of observations
  int n = x.size();
  
  // Get the number of knots and tau parameters
  int n_knots = sum(knots_mult);
  int n_knots_0 = knots_mult.size();
  int n_tau = (sum(knots_mult) - 1) - degree;
  
  // Assign tau parameters
  NumericVector tau = NumericVector(n_tau);
  for(int i = 0; i < n_tau; i++)
  {
    tau[i] = par[i];
  }
  
  // Assign knots
  NumericVector knots_0 = NumericVector(n_knots_0);
  knots_0[0] = min(x);
  knots_0[n_knots_0 - 1] = max(x);
  for(int i = 1; i < (n_knots_0 - 1); i++)
  {
    knots_0[i] = par[i - 1 + n_tau];
  }

  // Assign knots considering their multiplicity
  NumericVector knots = NumericVector(n_knots);
  int knots_counter = 0;
  for(int i = 0; i < n_knots; i++)
  {
    for(int j = 0; j < knots_mult[i]; j++)
    {
      knots[knots_counter] = knots_0[i];
      knots_counter++;
    }
  }
  
  // Validate the distance between the knots
  // and their sorting
  NumericVector knots_ecdf= NumericVector(n_knots_0);
  knots_ecdf[0] = 1 / n;
  knots_ecdf[n_knots_0 - 1] = 1;
  for(int i = 1; i < (n_knots_0 - 1); i++)
  {
    knots_ecdf[i] = mean(knots_0[i] >= x);
  }

  for(int i = 0; i < (n_knots_0 - 1); i++)
  {
    if((knots_ecdf[i + 1] - knots_ecdf[i]) < knots_dist)
    {
      return(1e+20);
    }
  }

  // Deal with b-splines
  List b = hpa::bsplineGenerate(knots, degree, false);
  List b_new = hpa::bsplineComb(b, tau);
  NumericMatrix b_m = b_new["m"];
  
  // Moments of SPL distribution
  double spl_exp = hpa::ehsa(b_m, knots, 0.0, 1.0, 1.0);
  double spl_exp2 = hpa::ehsa(b_m, knots, 0.0, 1.0, 2.0);
  double spl_sd = sqrt(spl_exp2 - std::pow(spl_exp, 2));
  
  // Estimate the density
  NumericVector den = hpa::dhsa(x, b_m, knots, 0.0, 1.0, false);
  den = den * spl_sd;

  // Minimize the distance
  NumericVector dist = den - den_kernel;
  dist = pow(dist, 2);
  double dist_value = sum(dist);

  return(dist_value);
}

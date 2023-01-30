# GARCH main function
npgarch <- function(y, x0 = NULL,
                    distribution = "normal",
                    arma_order = c(0, 0),
                    garch_order = c(1, 1),
                    degree = 1,
                    knots_mult = degree + 2,
                    gena_args = list(),
                    kernel_args = NULL)
{
  # Get initial x0 state
  is_x0_null <- is.null(x0)
  
  # Variable to store information on
  # preliminary estimation of knots
  kts <- NULL
  
  # Variable to store preliminary GARCH estimation
  garch_0 <- NULL
  
  # Some default values
  error_value <- -(1e+20)
  
  # Validation
  if (length(knots_mult) == 1)
  {
    knots_mult <- rep(1, knots_mult)
  }
  
  # Calculate mean and variance of 
  # the dependent variable
  y_mean <- mean(y)
  y_var <- var(y)
  
  # Set indexes of some variables
  ind <- list()
  n_par <- 0
  if (distribution != "KDE")
  {
    ind$mu <- 1
    ind$omega <- 2
    n_par <- 2
  }
  
  # Create vector to store parameters names
  par_names <- NULL 
  if (distribution != "KDE")
  {
    par_names <- c("mu", "omega")
  }
  
  # Indeсes of GARCH part
  if (garch_order[1] > 0)
  {
    ind$alpha <- (n_par + 1):(n_par + garch_order[1])
    n_par <- n_par + length(ind$alpha)
    par_names <- c(par_names, paste0("alpha", 1:garch_order[1]))
  }
  ind$beta <- NULL
  if (garch_order[2] > 0)
  {
    ind$beta <- (n_par + 1):(n_par + garch_order[2])
    n_par <- n_par + length(ind$beta)
    par_names <- c(par_names, paste0("beta", 1:garch_order[2]))
  }
  
  # Indexes of ARMA part
  ind$ar <- NULL
  if (arma_order[1] > 0)
  {
    ind$ar <- (n_par + 1):(n_par + arma_order[1])
    n_par <- n_par + length(ind$ar)
    par_names <- c(par_names, paste0("ar", 1:arma_order[1]))
  }
  ind$ma <- NULL
  if (arma_order[2] > 0)
  {
    ind$ma <- (n_par + 1):(n_par + arma_order[2])
    n_par <- n_par + length(ind$ma)
    par_names <- c(par_names, paste0("ma", 1:arma_order[2]))
  }

  # Assign initial values for optimization if need
  if (is.null(x0))
  {
    x0 <- rep(0, n_par + ifelse(distribution == "KDE", 2, 0))
    x0[ind$mu] <- y_mean
    x0[ind$omega] <- y_var
    if ((distribution == "SPL") | 
        (distribution == "PGN") |
        (distribution == "KDE"))
    {
      garch_0 <- npgarch::npgarch(y = y, x0 = x0,
                                  distribution = "normal",
                                  garch_order = garch_order,
                                  arma_order = arma_order,
                                  gena_args = gena_args)
      if (distribution != "KDE")
      {
        x0[c(ind$mu, ind$omega, 
             ind$alpha, ind$beta, 
             ind$ar, ind$ma)] <- garch_0$par
      }
      else
      {
        x0 <- rep(0, n_par)
        x0[c(ind$alpha, ind$beta, 
             ind$ar, ind$ma)] <- garch_0$par[-c(1, 2)]
      }
    }
  }

  if (distribution == "SPL")
  {
    # Get the number of parameters
    # related to spline
    n_knots <- length(knots_mult)
    n_tau <- (sum(knots_mult) - 1) - degree
    # Assign indeсes related to spline
    ind$tau <- (n_par + 1):(n_par + n_tau)
    n_par <- n_par + length(ind$tau)
    par_names <- c(par_names, paste0("tau", 1:length(ind$tau)))
    ind$inner_knots <- (n_par + 1):(n_par + n_knots - 2)
    n_par <- n_par + length(ind$inner_knots)
    par_names <- c(par_names, paste0("knot", 1:length(ind$inner_knots)))
    # Assign initial values to parameters 
    # related to spline
    if (is_x0_null)
    {
      ind_Rcpp <- list(mu = ind$mu, omega = ind$omega, 
                       alpha = ind$alpha, beta = ind$beta,
                       ar = ind$ar, ma = ind$ma)
      ind_Rcpp <- lapply(ind_Rcpp, function(x){x - 1})
      eps_tmp <- npgarch::GARCH_lnL(par = x0[c(ind$mu, ind$omega,
                                               ind$alpha, ind$beta,
                                               ind$ar, ind$ma)], 
                                    y = y, 
                                    distribution = "normal", 
                                    arma_order = arma_order,
                                    garch_order = garch_order,
                                    return_type = "eps",
                                    ind = ind_Rcpp)
      kts <- npgarch::KernelToSpline(x = eps_tmp$eps_std, 
                                     degree = degree, 
                                     knots_mult = knots_mult)
      x0 <- c(x0, rep(0, length(ind$tau) + length(ind$inner_knots)))
      x0[ind$tau] <- kts$tau
      x0[ind$inner_knots] <- kts$inner_knots
    }
  }
  if (distribution == "KDE")
  {
    ind$tau <- n_par + 1
    n_par <- n_par + 1
  }

  if (distribution == "PGN")
  {
    # Assign indeсes
    ind$tau <- (n_par + 1):(n_par + degree)
    n_par <- n_par + length(ind$tau)
    if (is_x0_null)
    {
      x0[ind$tau] <- 1e-16
    }
    par_names <- c(par_names, paste0("tau", 1:length(ind$tau)))
  }


  if (distribution == "KDE3")
  {
    if (is_x0_null)
    {
      bw0 <- density(garch_0$eps_std, bw = "SJ")$bw
      x0 <- c(x0, bw0)
    }
  }
  
  if (distribution == "KDE")
  {
    if (is_x0_null)
    {
      x0[ind$tau] <- density(garch_0$eps_std, bw = "SJ")$bw
    }
  }

  # Convert indeсes to Rcpp indexes
  ind_Rcpp <- lapply(ind, function(x){x - 1})

  # Genetic algorithm
  
    # Deal with initial arguments
  gena_args$fn <- npgarch::GARCH_aggregate
  #gena_args$constr.method <- npgarch::GARCH_constr
  
  gena_args$y <- y
  gena_args$degree <- degree
  gena_args$distribution <- distribution
  gena_args$knots_mult <- knots_mult
  gena_args$kernel_args <- kernel_args
  gena_args$garch_order <- garch_order
  gena_args$arma_order <- arma_order
  gena_args$ind <- ind_Rcpp
  
  gena_args$error_value <- error_value
  
  if ((distribution == "normal") |
      (distribution == "PGN"))
  {
    if (all(arma_order == 0))
    {
      gena_args$gr = npgarch::GARCH_lnL_grad
    }
  }
  
    # Set bounds
  lower <- rep(NA, n_par)
  upper <- rep(NA, n_par)
  delta_bound <- 1e-10
  lower[ind$mu] <- min(y_mean * 0.8, y_mean * 1.2)
  upper[ind$mu] <- max(y_mean * 0.8, y_mean * 1.2)

  lower[ind$omega] <- 0
  upper[ind$omega] <- y_var * 100

  if (garch_order[1] > 0)
  {
    lower[ind$alpha] <- delta_bound
    upper[ind$alpha] <- 1 - delta_bound
  }
  
  if (garch_order[2] > 0)
  {
    lower[ind$beta] <- delta_bound
    upper[ind$beta] <- 1 - delta_bound
  }
  
  if (arma_order[1] > 0)
  {
    lower[ind$ar] <- -2
    upper[ind$ar] <- 2
  }
  
  if (arma_order[2] > 0)
  {
    lower[ind$ma] <- -2
    upper[ind$ma] <- 2
  }
  
  if (distribution == "PGN")
  {
    lower[ind$tau] <- rep(-10, degree)
    upper[ind$tau] <- rep(10, degree)
  }
  
  if (distribution == "SPL")
  {
    lower[ind$tau] <- -10 * abs(x0[ind$tau])
    upper[ind$tau] <- 10 * abs(x0[ind$tau])
    
    lower[ind$inner_knots] <- -10 * abs(x0[ind$inner_knots])
    upper[ind$inner_knots] <- 10 * abs(x0[ind$inner_knots])
  }
  

  if (distribution == "KDE3")
  {
    lower <- c(lower, 1e-8)
    upper <- c(upper, x0[5] * 2)
  }
  

  if (distribution == "KDE")
  {
    lower[ind$tau] <- 0
    upper[ind$tau] <- 1
  }

  gena_args$lower <- lower
  gena_args$upper <- upper
  
  if (is.null(gena_args$pop.initial))
  {
    gena_args$pop.initial <- x0
  }
  else
  {
    gena_args$pop.initial <- rbind(gena_args$pop.initial, x0)
  }
  
    # Start the genetic algorithm
  gena_results <- do.call(what = gena, args = gena_args)

    # Get the estimates
  par_est <- gena_results$par
  names(par_est) <- par_names
  lnL_est <- gena_results$value
  
  # Estimate conditional volatility
  sigma <- npgarch::GARCH_lnL(par = par_est, 
                              y = y, degree = degree,
                              distribution = distribution,
                              arma_order = arma_order,
                              garch_order = garch_order,
                              knots_mult = knots_mult, 
                              return_type = "sigma",
                              error_value = error_value,
                              ind = ind_Rcpp)$val
  
  eps_list <- npgarch::GARCH_lnL(par = par_est, 
                                 y = y, degree = degree,
                                 distribution = distribution,
                                 arma_order = arma_order,
                                 garch_order = garch_order,
                                 knots_mult = knots_mult, 
                                 return_type = "eps",
                                 error_value = error_value,
                                 ind = ind_Rcpp)
  
  # Store various coefficients
  coef <- list()
  coef$mu <- par_est[ind$mu]
  coef$omega <- par_est[ind$omega]
  if (length(ind$alpha) > 0)
  {
    coef$alpha <- par_est[ind$alpha]
  }
  if (length(ind$beta) > 0)
  {
    coef$beta <- par_est[ind$beta]
  }
  if (length(ind$ar) > 0)
  {
    coef$ar <- par_est[ind$ar]
  }
  if (length(ind$ma) > 0)
  {
    coef$ma <- par_est[ind$ma]
  }
  if (length(ind$tau) > 0)
  {
    coef$tau <- par_est[ind$tau]
  }
  if (length(ind$inner_knots) > 0)
  {
    coef$inner_knots <- par_est[ind$inner_knots]
  }
  
  # Aggregate the results
  return_list <- list(par = par_est,
                      coef = coef,
                      y = y,
                      logLikelihood = lnL_est,
                      AIC = 2 * (n_par - lnL_est),
                      BIC = log(n) * n_par - 2 * lnL_est,
                      HQC = 2 * (log(log(n)) * n_par - lnL_est),
                      sigma = sigma,
                      gena = gena_results,
                      ind = ind,
                      ind_Rcpp = ind_Rcpp,
                      distribution = distribution,
                      knots_mult = knots_mult,
                      arma_order = arma_order,
                      garch_order = garch_order,
                      degree = degree,
                      eps = eps_list$eps,
                      eps_std = eps_list$eps_std)
  
  if (distribution == "PGN")
  {
    return_list$tau <- par_est[ind$tau]
  }
  
  if (distribution == "SPL")
  {
    return_list$tau <- par_est[ind$tau]
    return_list$inner_knots <- par_est[ind$inner_knots]
    return_list$knots <- c(min(eps_list$eps_std), 
                           return_list$inner_knots,
                           max(eps_list$eps_std))
  }
  
  class(return_list) <- "npgarch"
  return(return_list)                             
}

print.npgarch <- function(object, ...)
{
  # Get some variables
  ind <- object$ind
  ind_Rcpp <- object$ind_Rcpp
  y <- object$y
  AIC <- object$AIC
  BIC <- object$BIC
  HQC <- object$HQC
  sigma <- object$sigma
  par <- object$par
  arma_order <- object$arma_order
  garch_order <- object$garch_order
  distribution <- object$distribution
  degree <- object$degree
  knots_mult <- object$knots_mult
  
  # Estimate covariance matrix
  H <- gena::gena.hessian(fn = GARCH_aggregate, par = par,
                          fn.args = list(distribution = distribution,
                          knots_mult = knots_mult,
                          degree = degree,
                          y = y, ind = ind_Rcpp, 
                          arma_order = arma_order,
                          garch_order = garch_order))
  cov_mat <- qr.solve(-H, tol = 1e-16)
  colnames(cov_mat) <- names(par)
  rownames(cov_mat) <- names(par)
  se <- sqrt(diag(cov_mat))
  
  # Calculate p-values
  z_value <- par / se
  p_value <- rep(NA, length(par))
  for (i in 1:length(par))
  {
    p_value[i] <- 2 * min(pnorm(z_value[i]), 1 - pnorm(z_value[i]))
  }
  tbl <- cbind(par, se, p_value)
  
  return(tbl)
}

GARCH_constr <- function(population,
                         method = "no",
                         par,
                         iter,
                         ...)
{

  if (distribution == "KDE")
  {
    delta <- 0.001
    
    par[1] <- min(max(par[1], delta), 1 - delta)
    par[2] <- min(max(par[2], delta), 1 - delta)
    
    if (par[3] <= 0)
    {
      par[3] <- delta
    }
    
    return(par)
  }
  
  delta <- 0.001
  
  par[ind$omega] <- max(0, par[2])
  par[ind$alpha] <- min(max(par[ind$alpha], delta), 1 - delta)
  par[ind$beta] <- min(max(par[ind$beta], delta), 1 - delta)
  
  if ((par[ind$alpha] + par[ind$beta]) >= 1)
  {
    coef_sum <- par[ind$alpha] + par[ind$beta]
    
    par[ind$alpha] <- (par[ind$alpha] / coef_sum) - delta
    par[ind$beta] <- (par[ind$beta] / coef_sum) - delta
  }

  if (distribution == "SPL")
  {
    n_knots <- sum(knots_mult)
    n_tau <- (n_knots - 1) - degree
    
    par[ind$inner_knots] <- sort(par[ind$inner_knots])
  }
  

  if (distribution == "KDE3")
  {
    if (par[5] <= 0)
    {
      par[5] <- 0.0000001
    }
  }

  return(par)
}
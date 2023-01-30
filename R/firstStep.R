#
KernelToPolynomial <- function(x,
                               K,
                               den_est)
{
  ker <- approxfun(density(x, bw = "SJ"))
  den_est <- ker(x)
  
  opt <- optim(par = rep(0, K), 
               fn = KernelToPolynomial_min, 
               method = "Nelder-Mead",
               control = list(maxit = 1000000000,
                              reltol = 1e-10),
               x = x, 
               den_est = den_est)
  
  return(opt$par)
}

#
KernelToSpline <- function(x,
                           degree,
                           knots_mult)
{
  n <- length(x)
  
  n_knots <- length(knots_mult)
  n_tau <- (sum(knots_mult) - 1) - degree
  tau <- rep(1, n_tau)
  knots_0 <- as.numeric(quantile(x, seq(0, 1, length = n_knots)))
  
  ker <- approxfun(density(x, bw = "SJ"))
  den_kernel <- ker(x)
  
  x0 <- c(tau, knots_0[-c(1, n_knots)])
  opt <- optim(par = x0, 
               fn = npgarch::KernelToSpline_LS, 
               method = "Nelder-Mead",
               control = list(maxit = 1000000000,
                              reltol = 1e-10,
                              abstol = 1e-10),
               x = x, 
               degree = degree, 
               knots_mult = knots_mult,
               den_kernel = den_kernel)

  inner_knots <- opt$par[-(1:n_tau)]
  boundary_knots <- c(min(x), max(x))
  knots <- c(boundary_knots[1], inner_knots, boundary_knots[2])
  tau <- opt$par[1:n_tau]

  val <- list(knots = knots,
              inner_knots = inner_knots,
              boundary_knots = boundary_knots,
              tau = tau,
              value = opt$val / n)
  
  return(val)
}
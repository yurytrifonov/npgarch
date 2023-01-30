# Function to simulate GARCH process
# given appropriate standardized
# shocks' distribution
GARCH.Simulate <- function (n,                         # number of observations
                            mu, omega, alpha, beta,    # parameters
                            rden,                      # function to simulate
                                                       # standardized shocks   
                            rden_args)                 # list of rden arguments                      
{
  # simulate standardized shocks
  rden_args$n <- n
  xi <- do.call(rden, rden_args)
  
  # initialize main variables
  eps <- rep(NA, n)                               # shocks
  sigma <- rep(NA, n)                             # conditional volatility
  
  # set initial values
  sigma_0 <- sqrt(omega / (1 - alpha - beta))     # unconditional variance
  sigma[1] <- sigma_0                             # conditional volatility at 
                                                  # the first period
  eps[1] <- sigma[1] * xi[1]                      # shock in the first period
  
  # calculate conditional volatility
  for(t in 2:n)
  {
    sigma[t] <- sqrt(omega +
                     alpha * eps[t - 1] ^ 2 + 
                     beta * sigma[t - 1] ^ 2)
    eps[t] <- sigma[t] * xi[t]
  }
  
  # aggregate the results
  df <- data.frame(sigma = sigma, 
                   y = eps + mu,
                   eps = eps,
                   eps_std = eps / sigma)
  
  return(df)
}

#
spl.crossover <- function(children,
                          n_pop,
                          n_genes,
                          n_tau,
                          n_knots,
                          boundary_knots,
                          prob = 0.8)
{
  ind_1 <- 4 + n_tau                         # for mean, omega, alpha, beta and tau parameters
  ind_2 <- (n1 + 2):(n_genes - 1)            # for knots parameters
  
  for (i in 1:n_pop)
  {
    ind <- c(i - 1, i)
    
    # first
    weight <- runif(n = n_genes, min = 0, max = 1)
    children[ind[1], ind_1] <- parents[ind[1], ind_1] * weight[ind_1] +
                               parents[ind[2], ind_1] * (1 - weight[ind_1])
    children[ind[2], ind_1] <- parents[ind[2], ind_1] * weight[ind_1] +
                               parents[ind[1], ind_1] * (1 - weight[ind_1])
    
    # second
    knots_diff_1 <- diff(c(boundary_knots[ind[1], 1], 
                           parents[ind[1], ind_2]))
    knots_diff_2 <- diff(c(boundary_knots[ind[2], 1], 
                           parents[ind[2], ind_2]))
    
    knots_diff_c1 <- knots_diff_1 * weight[ind_2] +
                     knots_diff_2 * (1 - weight[ind_2])
    knots_diff_c2 <- knots_diff_2 * weight[ind_2] +
                     knots_diff_1 * (1 - weight[ind_2])
    
    
    knots_c1 <- sample(1, c(boundary_knots[ind, 1])) + cumsum(knots_diff_c1)
    knots_c2 <- sample(1, c(boundary_knots[ind, 1])) + cumsum(knots_diff_c2)
    
    children[ind[1], ind_2] <- knots_c1
    children[ind[2], ind_2] <- knots_c2
  }
  
  return(children)
}

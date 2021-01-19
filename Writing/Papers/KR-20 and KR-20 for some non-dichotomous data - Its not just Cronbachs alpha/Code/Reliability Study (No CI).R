library(MCMCpack) # Necessary for rinvgamma()
library(psych) 
library(MBESS)

sim.data <- function(rho, n_items, n_subjects, mu, distribution) {
  
  # This function simulates normal and non-normal test data with a population test
  # reliability given by rho. 
  
  # Inputs: 
  # rho = Population test reliability
  # n_items = Number of test items, the test length
  # n_subjects 
  # mu = Mean of underlying ability distribution. This does not directly affect the 
  #      population reliability rho but may affect estimation of rho by making the
  #      generating distributions more or less similar to a normal distribution.
  
  
  # Find the value of M necessary for population value rho given a test of length n_items.
  # In general, the relation is rho = n_items/(M + n_items)
  
  M = (1-rho)/rho*n_items
  
  # Create matrix of test item responses y
  
  y <- matrix(0, n_subjects, n_items)
  
  # Simulate data
  
  if(distribution == 'normal' | distribution == 'Normal') {
    
    sigma = 1
    
    theta = rnorm(n_subjects, mu, sqrt(1/M))
    
    for(i in 1:n_subjects) y[i,] = rnorm(n_items, theta[i], 1)
    
  }
  else if(distribution == 'binomial' | distribution == 'Binomial') {
    
    alpha = mu*M
    beta = (1-mu)*M
    
    if(mu < 0 | mu > 1) {
      print("Error: mu must be between 0 and 1 for the binomial distribution")
      return(NA)
    }
    
    theta = rbeta(n_subjects, alpha, beta)
    
    for(i in 1:n_subjects) y[i,] = rbinom(n_items, 1, theta[i])
    
  }
  else if(distribution == 'poisson' | distribution == 'Poisson') {
    
    alpha = mu*M
    beta = M
    
    if(alpha < 0 | beta < 0) {
      print("Error: Impossible values for mu and M, mu must be positive")
      return(NA)
    }
    
    theta = rgamma(n_subjects, alpha, beta)
    
    for(i in 1:n_subjects) y[i,] = rpois(n_items, theta[i])
    
  }
  else if(distribution == 'gamma' | distribution == 'Gamma') {
    
    alpha = M + 1
    beta = mu*M
    
    if(alpha < 0 | beta < 0) {
      print("Error: Impossible values for mu and M, mu must be positive")
      return(NA)
    }
    
    theta = rinvgamma(n_subjects, alpha, beta)
    
    for(i in 1:n_subjects) y[i,] = rexp(n_items, 1/theta[i])
    
    
  }
  else if(distribution == 'negative binomial' | distribution == 'Negative Binomial') {
    
    # Note that not all mu values are possible for a given M.
    # When M is large, mu should be chosen close to (but larger than) 1. 
    
    d1 = -2*mu*(M-1)/((mu*(mu-1)*M-2))
    d2 = 2*mu/(mu-1)
    
    # Check to ensure d1 > 0 and d2 > 4. If not, decrease mu to be closer to (but larger than) 1.
    
    if(d1 < 0 | d2 < 4) {
      print("Error: Impossible values for mu and M, move mu closer to (but larger than) 1")
      return(NA)
    }
    
    theta = rf(n_subjects, d1, d2)
    
    for(i in 1:n_subjects) y[i,] = rgeom(n_items, 1/(1 + theta[i])) # p = 1/(1 + theta)
    
    
  }
  else{
    print("Error, input one of the following distributions: normal, binomial, poisson, gamma, negative binomial")
    return(NA)
  }
  
  return(y)
  
}



#cronbach.alpha = function(y) dim(y)[2]/(dim(y)[2]-1)*(1 - sum(apply(y,2,var))/var(apply(y,1,sum)))

get.reliability.estimates <- function(y, distribution) {
  
  # Input: y, matrix of test responses. Should be n_subjects rows and n_items columns.
  
  n_subjects = dim(y)[1]
  n_items = dim(y)[2]
  
  x = apply(y,1,sum)
  
  # Cronbach's alpha
  
  alpha_est = dim(y)[2]/(dim(y)[2]-1)*(1 - sum(apply(y,2,var))/var(apply(y,1,sum)))
  
  # Omega
  
  omega_est = omega(y, plot = F)
  
  
  # KR20
  
  kr20 = NA
  
  if(distribution == 'binomial' | distribution == 'Binomial') {kr20 = n_items/(n_items - 1)*(1 - (sum(apply(y,2,mean)*(1-apply(y,2,mean)))/var(x)))}
  if(distribution == 'poisson' | distribution == 'Poisson') {kr20 = (1 - (sum(apply(y,2,mean))/var(x)))}
  if(distribution == 'gamma' | distribution == 'Gamma') {kr20 = n_items/(n_items + 1)*(1 - (sum(apply(y,2,mean)^2)/var(x)))}
  if(distribution == 'negative binomial' | distribution == 'Negative Binomial') {kr20 = n_items/(n_items + 1)*(1 - (sum(apply(y,2,mean) + apply(y,2,mean)^2)/var(x)))}
  
  
  # KR21
  
  kr21 = NA
  
  if(distribution == 'binomial' | distribution == 'Binomial') {kr21 = n_items/(n_items - 1)*(1 - n_items*mean(x/n_items)*(1 - mean(x/n_items))/var(x))}
  if(distribution == 'poisson' | distribution == 'Poisson') {kr21 = (1 - mean(x)/var(x))}
  if(distribution == 'gamma' | distribution == 'Gamma') {kr21 = n_items/(n_items + 1)*(1 - n_items*mean(x/n_items)^2/var(x))}
  if(distribution == 'negative binomial' | distribution == 'Negative Binomial') {kr21 = n_items/(n_items + 1)*(1 - n_items*(mean(x/n_items) + mean(x/n_items)^2)/var(x))}
  
  # GLB
  
  glb.alg = glb.algebraic(cov(y))$glb
  
  return(list('alpha' = alpha_est,'omega.h' = omega_est$omega_h, 'omega.tot' = omega_est$omega.tot, 'kr21' = kr21, 'kr20' = kr20, 'glb' = glb.alg))
  
  
}


simulation.study <- function(n_sims, rho, n_items, n_subjects, mu, distribution) {
  
  alpha = rep(0, n_sims)
  omega.h = rep(0, n_sims)
  omega.tot = rep(0, n_sims)
  kr20 = rep(0, n_sims)
  kr21 = rep(0, n_sims)
  glb = rep(0, n_sims)
  
  start.time = Sys.time()
  
  for(sim in 1:n_sims) {
    
    y <- sim.data(rho, n_items, n_subjects, mu, distribution)
    
    est <- get.reliability.estimates(y, distribution)
    
    alpha[sim] = est$alpha
    omega.h[sim] = est$omega.h
    omega.tot[sim] = est$omega.tot
    kr20[sim] = est$kr20
    kr21[sim] = est$kr21
    glb[sim] = est$glb
    
    if(sim %% 1e5 == 0) {
      curr.time = Sys.time()
      print(paste(sim, "/", n_sims, ", Time Taken = ", curr.time - start.time, sep=""))
    }
    
    
  }
  
  return(list('n_sims' = n_sims, 'n_subjects' = n_subjects, 'n_items' = n_items, 'rho' = rho, 'mu' = mu, 'distribution' = distribution, 'alpha' = alpha, 'omega.h' = omega.h, 'omega.tot' = omega.tot, 'kr20' = kr20, 'kr21' = kr21, 'glb' = glb))
  
}


# Simulation Study


n_sims = 1e3
rho = 0.8
n_items = 5
n_subjects = 100
mu = 0.5
distribution = 'binomial'



start = Sys.time()

res <- simulation.study(n_sims, rho, n_items, n_subjects, mu, distribution)

stop = Sys.time()

stop - start

alpha.rmse = round(sqrt(mean((res$alpha - rho)^2)),3)
kr20.rmse = round(sqrt(mean((res$kr20 - rho)^2)),3)
kr21.rmse = round(sqrt(mean((res$kr21 - rho)^2)),3)
omega.h.rmse = round(sqrt(mean((res$omega.h - rho)^2)),3)
omega.tot.rmse = round(sqrt(mean((res$omega.tot - rho)^2)),3)
glb.rmse = round(sqrt(mean((res$glb - rho)^2)),3)

alpha.bias = round(mean(res$alpha - rho),3)
kr20.bias = round(mean(res$kr20 - rho),3)
kr21.bias = round(mean(res$kr21 - rho),3)
omega.h.bias = round(mean(res$omega.h - rho),3)
omega.tot.bias = round(mean(res$omega.tot- rho),3)
glb.bias = round(mean(res$glb - rho),3)

alpha.sd = round(sd(res$alpha),3)
kr20.sd = round(sd(res$kr20),3)
kr21.sd = round(sd(res$kr21),3)
omega.h.sd = round(sd(res$omega.h),3)
omega.tot.sd = round(sd(res$omega.tot),3)
glb.sd = round(sd(res$glb),3)

results.table = data.frame( 'Omega.H' = c(omega.h.rmse, omega.h.bias, omega.h.sd), 'Omega.tot' = c(omega.tot.rmse, omega.tot.bias, omega.tot.sd), 'alpha' = c(alpha.rmse, alpha.bias, alpha.sd), 'kr20' = c(kr20.rmse, kr20.bias, kr20.sd),  'kr21' = c(kr21.rmse, kr21.bias, kr21.sd), 'glb' = c(glb.rmse, glb.bias, glb.sd))
row.names(results.table) = c('RMSE', 'Bias', 'SD')
results.table

stop - start

# Cronbach's alpha and KR20/21 simulation
# code for the Binomial - Beta model

# Specify simulation parameters 

n_items = 5
n_subjects = 30

true_mu = .5

true_rho = 0.8 # This is the true population reliability

n_sims = 1e5

# Calculate the M required to achieve the population reliability
# given the test length. The formula is 
# true_rho = n_items/(true_M + n_items)

true_M = (1-true_rho)/true_rho*n_items

# Set the seed for consistent results

set.seed(5)

# Functions

cronbach.alpha = function(y) {

  # Calculates Cronbach's alpha for an input matrix of scores y
  # where each row gives the item scores for a single subject

  n = dim(y)[2]
  return(n/(n-1)*(1 - sum(apply(y,2,var))/var(apply(y,1,sum))))
}


# Create vectors to store output

alpha = rep(0, n_sims)
kr.20 = rep(0, n_sims)
kr.21 = rep(0, n_sims)

# Begin loop over all simulations

for(sim in 1:n_sims) {
  
  # Simulate data
  
  # Generate true abilities

  # For beta distribution the following transformations are used:
  #
  # mu = alpha/(alpha + beta)
  # M = alpha + beta
  #
  # alpha = M*mu
  # beta = M*(1-mu)
  
  theta <- rbeta(n_subjects, true_M*true_mu, true_M*(1-true_mu))
  
  # Create scored response matrix and simulate scores for each subject using a binomial distribution
  
  y <- matrix(0, n_subjects, n_items)
  
  for(i in 1:n_subjects) y[i,] = rbinom(n_items, 1, theta[i])
  
  # Get raw sum scores for each subject
  
  x = apply(y,1,sum)
  
  # Calculate the generalized kr-20 and kr-21 formulas
  
  kr.20[sim] = n_items/(n_items - 1)*(1 - (sum(apply(y,2,mean)*(1-apply(y,2,mean)))/var(x)))
  kr.21[sim] = n_items/(n_items - 1)*(1 - n_items*mean(x/n_items)*(1 - mean(x/n_items))/var(x))
  
  # Cronbach's alpha

  alpha[sim] = cronbach.alpha(y)

  # A counter to monitor progress, if desired
  
  #if(sim %% 1e5 == 0) print(paste('Sim ', sim, '/', n_sims, ' Complete', sep = ''))
  
}

# Calculate RMSE, bias, and standard deviation for each estimator

alpha.rmse = round(sqrt(mean((alpha - true_rho)^2)),3)
alpha.bias = round(mean(alpha - true_rho),3)
alpha.sd = round(sd(alpha),3)

kr.20.rmse = round(sqrt(mean((kr.20 - true_rho)^2)),3)
kr.20.bias = round(mean(kr.20 - true_rho),3)
kr.20.sd = round(sd(kr.20),3)

kr.21.rmse = round(sqrt(mean((kr.21 - true_rho)^2)),3)
kr.21.bias = round(mean(kr.21 - true_rho),3)
kr.21.sd = round(sd(kr.21),3)

# Put the results into a table and print it

res = data.frame('alpha' = c(alpha.rmse, alpha.bias, alpha.sd), 'kr20' = c(kr.20.rmse, kr.20.bias, kr.20.sd), 'kr21' = c(kr.21.rmse, kr.21.bias, kr.21.sd))
row.names(res) = c('RMSE', 'Bias', 'SD')
res




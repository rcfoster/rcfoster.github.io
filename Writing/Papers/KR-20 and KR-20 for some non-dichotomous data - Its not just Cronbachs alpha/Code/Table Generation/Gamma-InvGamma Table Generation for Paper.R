# This code is to generate Table 3, simulation results for the gamma - inverse gamma model, in the paper. 
# It is essentially the file "Gamma-InvGamma generalized KR simulation.R" but without the 
# additional comments and wrapped in a loop over all possible simulation values. All of the
# data is thrown away after each loop and only LaTeX code for each row of the the interior of 
# the table is preserved.


library(MCMCpack) # Necessary for rinvgamma()


# Functions

cronbach.alpha = function(y) {
  n = dim(y)[2]
  return(n/(n-1)*(1 - sum(apply(y,2,var))/var(apply(y,1,sum))))
}


# Specify simulation parameters and generate list of parameters for simulation


n_items_list = c(5, 10, 30)
n_subjects_list = c(30, 75, 500)

true_mu = 1

true_rho_list = c(0.3, 0.6, 0.8)

pars <- NULL

for(i in 1:length(true_rho_list)) {
  for(j in 1:length(n_items_list)) {
    for(k in 1:length(n_subjects_list)) {
      pars <- rbind(pars, data.frame('mu' = true_mu, 'true_rho' = true_rho_list[i], 'n_items' = n_items_list[j], 'n_subjects' = n_subjects_list[k], 'true_M' = (1-true_rho_list[i])/true_rho_list[i]*n_items_list[j]))
      
    }
  }
}


n_sims = 1e6

tab <- NULL

start = Sys.time()

for(par in 1:dim(pars)[1]) {
  
  # Set the seed for consistent results
  
  set.seed(5)
  
  # Create vectors to store output
  
  alpha = rep(0, n_sims)
  kr.20 = rep(0,n_sims)
  kr.21 = rep(0,n_sims)
  
  # Get parameters
  
  true_mu = pars[par, 1]
  true_rho = pars[par, 2]
  n_items = pars[par,3]
  n_subjects = pars[par,4]
  true_M = pars[par,5]
  
  # Begin loop over all simulations
  
  for(sim in 1:n_sims) {
    
    # Simulate data
    
    #alpha = true_M + 1
    #beta = mu*M
    
    theta <- rinvgamma(n_subjects, true_M+1, true_mu*true_M)
    
    # Create scored response matrix 
    
    y <- matrix(0, n_subjects, n_items)
    
    for(i in 1:n_subjects) y[i,] = rexp(n_items, 1/theta[i])
    
    # Get the vector of response lengths and raw scores for each subject
    
    x = apply(y,1,sum)
    n = rep(n_items, n_subjects)
    
    # Calculate the generalized kr-20 and kr-21 formulas
    
    kr.20[sim] = n_items/(n_items + 1)*(1 - (sum(apply(y,2,mean)^2)/var(x)))
    kr.21[sim] = n_items/(n_items + 1)*(1 - n_items*mean(x/n_items)^2/var(x))
    
    # Cronbach's alpha
    
    alpha[sim] = cronbach.alpha(y)
     
  }
  
  alpha.rmse = format(round(sqrt(mean((alpha - true_rho)^2)),3),nsmall = 3)
  alpha.bias = format(round(mean(alpha - true_rho),3), nsmall = 3)
  alpha.sd = format(round(sd(alpha),3), nsmall = 3)
  
  kr.20.rmse = format(round(sqrt(mean((kr.20 - true_rho)^2)),3),nsmall = 3)
  kr.20.bias = format(round(mean(kr.20 - true_rho),3),nsmall = 3)
  kr.20.sd = format(round(sd(kr.20),3),nsmall = 3)
  
  kr.21.rmse = format(round(sqrt(mean((kr.21 - true_rho)^2)),3),nsmall = 3)
  kr.21.bias = format(round(mean(kr.21 - true_rho),3),nsmall = 3)
  kr.21.sd = format(round(sd(kr.21),3),nsmall = 3)
  
  table.row = paste(true_rho, " & ", n_subjects, " & ", n_items, " & ", true_mu, " & ", round(true_M,2), " & ", alpha.rmse, " (", alpha.bias, ", ", alpha.sd, ") & ", kr.20.rmse, " (", kr.20.bias, ", ", kr.20.sd, ") & ", kr.21.rmse, " (", kr.21.bias, ", ", kr.21.sd, ") \\", sep = '')
  
  tab <- rbind(tab, table.row)

  curr = Sys.time()
  
  print(paste("Simulation ", par, "/", dim(pars)[1], " complete, Time taken = ", curr - start, sep = ''))
  
  
}

rownames(tab) = rep('', dim(tab)[1])
print(tab, quote = F)
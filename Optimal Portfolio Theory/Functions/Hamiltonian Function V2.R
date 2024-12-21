

###########Generalised_network_code##############

hamil_func <- function(N, S, lambda_0, n_obs, Omega_init, tau_init, epsilons, A, catrate, mu_beta, var_beta, step_size, n_steps){
  # lambda1 ~ G(a_0, b_0)
  p <- ncol(Omega_init)
  Omega_chain <- array(NA, dim = c(N+1, p, p))
  tau_chain <- array(NA, dim = c(N+1, p, p))
  lambda_1_chain <- array(NA, dim = c(N+1, p, p))
  beta_chain <- array(NA, dim = c(length(A) + 1,N+1))## A is a list - this is the numbe rof list items
  Omega_chain[1,,] <- Omega_init
  Sigma <- pd.solve(Omega_chain[1,,])
  tau_chain[1,,] <- tau_init
  #NEW
  for(i in 1:length(epsilons)){
    beta_chain[i, 1] <- mu_beta[i]
  }
  accrate <- numeric(N) #we are now updating both betas at the same time, so no need for N X 2 dimensions, just N. 
  lambda_1_chain[1,,] <- exp(beta_chain[1,1] + Reduce(`+`, lapply(seq_along(A), function(j) beta_chain[j+1, 1] * A[[j]])))
  
  for(t in 1:N){
    
    #print('new iteration!!!')
    
    Omega_chain[t+1,,] <- Omega_chain[t,,]
    tau_chain[t+1,,] <- tau_chain[t,,]
    lambda_1_chain[t+1,,] <- lambda_1_chain[t,,]
    
    for(i in 1:length(epsilons)){
      beta_chain[i, t+1] <- beta_chain[i, t]
    }
    
    #NEW
    for(j in 1:p){
      v <- rgamma(1, shape = n_obs/2 + 1, rate = (S[j,j]+lambda_0)/2)
      # omega_j | tau, lambda
      #inv_Omega_mjmj <- pd.solve(Omega_chain[t+1,-j,-j])
      inv_Omega_mjmj <- Sigma[-j, -j] - Sigma[-j, j] %*% t( Sigma[-j, j] ) / Sigma[j, j]
      #C <- pd.solve((S[j,j]+lambda_chain[t+1])*inv_Omega_mjmj + diag(1/tau_chain[t+1, j, -j], p-1))## Can we do this quickly?! ## Woodbury
      C_inv <- (S[j,j]+lambda_0)*inv_Omega_mjmj + diag(1/tau_chain[t+1, j, -j], p-1) 
      chol_C_inv <- chol(C_inv)
      C <- chol2inv(chol_C_inv)
      omega_j <- rmvnorm(1, -C%*%S[j,-j], C)
      omega_j <- matrix(omega_j, nrow = 1)
      Omega_chain[t+1,j,-j] <- omega_j
      Omega_chain[t+1,-j,j] <- omega_j
      Omega_chain[t+1, j, j] <- v + omega_j%*%inv_Omega_mjmj%*%t(omega_j)
      # tau | omega_j, lambda
      tau_j <- 1/rinvgauss(p-1, sqrt((lambda_1_chain[t+1, j, -j]^2)/(Omega_chain[t+1, j, -j]^2)), lambda_1_chain[t+1, j, -j]^2 )
      tau_chain[t+1, j,-j] <- tau_j
      tau_chain[t+1, -j,j] <- tau_j
      
      Sigma <- update_Sigma(Sigma, i=j, Omega_11_inv = inv_Omega_mjmj, Omega_11_inv_X_omega_12 = inv_Omega_mjmj%*%t(omega_j), gam = drop(v))## does this need a negative or not?
    }
    
    print(beta_chain[,t+1])
    
    iterations <- 100
    
    beta_subchain <- array(NA, dim = c(iterations + 1,2))
    accrate_sub <- numeric(iterations)
    
    beta_subchain[1, ] <- beta_chain[, t+1]
    
    for (iter in 1:iterations){
    update_result <- hmc_update(beta_subchain[iter,], Omega_chain[t+1,,], A, mu_beta, var_beta, step_size = step_size, n_steps = n_steps)
    beta_subchain[iter + 1,] <- update_result$beta
    #print(beta_subchain[iter + 1,])
    accrate[iter] <- ifelse(update_result$accepted, 1, 0)
    #print(accrate[t])
    }
    
    beta_chain[, t+1] <- beta_subchain[iterations + 1,]
    
    lambda_1_chain[t+1,,] <- exp(beta_chain[1,t+1] + Reduce(`+`, lapply(seq_along(A), function(j) beta_chain[j+1, t+1] * A[[j]])))
    
    if((t %% (N/catrate)) == 10){
      cat("Iteration", t, "done", "\n")
    }
  }

  
  return(list("Omega_chain" =  Omega_chain, "lambda_1_chain" =  lambda_1_chain, 'beta_chain' = beta_chain, 'accrate' = accrate))
}







#####HMC update###########

hmc_update <- function(beta_vector, Omega, A, mu_beta, var_beta, step_size, n_steps) {
  current_energy <- -log_posterior(beta_vector, Omega, A, mu_beta, var_beta)
  current_momentum <- rnorm(length(beta_vector))  # Draw random momentum
  kinetic_energy <- sum(current_momentum^2) / 2
  
  # Simulate dynamics
  new_state <- leapfrog(beta_vector, current_momentum, Omega, A, mu_beta, var_beta, step_size, n_steps)
  new_beta_vector <- new_state$beta
  new_energy <- -log_posterior(new_beta_vector, Omega, A, mu_beta, var_beta)
  new_kinetic <- sum(new_state$momentum^2) / 2
  
  # Accept or reject the new state
  if (runif(1) < exp(current_energy + kinetic_energy - new_energy - new_kinetic)) {
    return(list(beta = new_beta_vector, accepted = TRUE))  # Accept
  } else {
    return(list(beta = beta_vector, accepted = FALSE))  # Reject
  }
}
#####HMC update###########


#######function for log likelihood computation#############

log_posterior <- function(beta_vector, Omega, A, mu_beta, var_beta) {
  # Assuming A is the adjacency matrix and Omega is the precision matrix
  #print(beta_chain[1,2])
  
  lambda_terms <- exp(beta_vector[1] + Reduce(`+`, lapply(seq_along(A), function(j) beta_vector[j+1] * A[[j]][lower.tri(A[[j]])])))
  laplace_likelihood <- sum(dlaplace(Omega[lower.tri(Omega, diag = FALSE)], 0, lambda_terms))
  normal_prior <- sum(dnorm(beta_vector[1], mu_beta, sqrt(var_beta), log = TRUE)) + sum(dnorm(beta_vector[2], mu_beta, sqrt(var_beta), log = TRUE))

  
  return(laplace_likelihood + normal_prior)
}


#######function for log likelihood computation#############




#######function for computation of log posterior gradient##############
gradient_log_posterior <- function(beta_vector, Omega, A, mu_beta, var_beta, epsilon = 1e-3) {
  grad <- numeric(length(beta_vector))
  
  for (i in seq_along(beta_vector)) {
    beta_eps <- beta_vector
    beta_eps[i] <- beta_eps[i] + epsilon
    grad[i] <- (log_posterior(beta_eps, Omega, A, mu_beta, var_beta) - log_posterior(beta_vector, Omega, A, mu_beta, var_beta)) / epsilon
    #print(paste('hopefully just a number:', grad[i]))
  }
  return(grad)
}

#######function for computation of log posterior gradient##############





#######function for leapfroging##############

leapfrog <- function(beta_vector, momentum, Omega, A, mu_beta, var_beta, step_size, n_steps) {

  # Make half step for momentum at the beginning
  momentum <- momentum - (step_size / 2) * gradient_log_posterior(beta_vector, Omega, A, mu_beta, var_beta)
  
  beta_proposal <- beta_vector
  
  for (step in 1:n_steps) {

    beta_proposal <- beta_proposal + step_size * momentum
    
    #print(step)
    
    # Make a full step for the momentum, except at end of trajectory
    if (step != n_steps) {
      momentum <- momentum - step_size * gradient_log_posterior(beta_proposal, Omega, A, mu_beta, var_beta)
    }
  }

  
  # Make half step for momentum at the end
  momentum <- momentum - (step_size / 2) * gradient_log_posterior(beta_proposal, Omega, A, mu_beta, var_beta)
  
  #print(paste('momentum', momentum))
  #print(paste('beta', beta_proposal))
  
  return(list(beta = beta_proposal, momentum = momentum))

  
}

#######function for leapfroging##############








# https://github.com/cran/ssgraph/blob/master/src/matrix.cpp
update_Sigma <- function(Sigma, i, Omega_11_inv, Omega_11_inv_X_omega_12, gam){#, p){
  #dim <- p
  #p1 <- dim - 1
  
  #alpha_gam = 1/gam
  #alpha_ij = - alpha_gam
  
  Sigma[ -i, -i ] = Omega_11_inv + Omega_11_inv_X_omega_12 %*% t( Omega_11_inv_X_omega_12 ) / gam
  
  Sigma[-i,i] = - Omega_11_inv_X_omega_12 / gam
  Sigma[i,-i] = - Omega_11_inv_X_omega_12 / gam
  
  Sigma[ i, i ]  = 1 / gam
  
  return(Sigma)
}



















Bayesian_GLASSO_prior_simulation_twoLambda <- function(p, lambda_1, lambda_0){
  # p is the dimension of the matrix 
  # lamda is the GLASSO parameter
  
  theta <- matrix(0, nrow = p, ncol = p)
  theta[lower.tri(theta, diag = FALSE)] <- rlaplace(p*(p-1)/2, 0, 1/lambda_1)
  theta <- theta + t(theta)
  diag(theta) <- rexp(p, lambda_0/2)
  return(theta)
}

Bayesian_GLASSO_hyperprior_simulation_twoLambda_lambda0Fixed <- function(p, a_1_0, b_1_0, lambda_0 = 2){
  # p is the dimension of the matrix 
  lambda_1 <- rgamma(1, shape = a_1_0, rate = b_1_0)
  theta <- Bayesian_GLASSO_prior_simulation_twoLambda(p, lambda_1, lambda_0)
  return(theta)
}

threshold <- function(Theta_mat, threshold){
  return(Theta_mat*(abs(Theta_mat) >= threshold))
}

Bayesian_GLASSO_hyperprior_summaries_lambda0Fixed <- function(p, a_1_0, b_1_0, lambda_0 = 2, N_MC = 1000, thresh = 0.01){
  pd_indicator <- rep(NA, N_MC)
  sparsity <- rep(NA, N_MC)
  for(j in 1:N_MC){
    theta_samp <- Bayesian_GLASSO_hyperprior_simulation_twoLambda_lambda0Fixed(p, a_1_0, b_1_0, lambda_0)
    pd_indicator[j] <- min(eigen(theta_samp)$values) >= 0
    rho_samp <- threshold(cov2cor(theta_samp), thresh)
    sparsity[j] <- sum(rho_samp == 0)/(p*(p-1))
  }
  return(list("pd_indicator" = mean(pd_indicator), "sparsity" = mean(sparsity))) 
}

posterior_sample_to_credibility_interval <- function(post_sample, alpha = 0.95){
  interval_lower <- quantile(post_sample, (1-alpha)/2)
  interval_upper <- quantile(post_sample, 1 - (1-alpha)/2)
  return(list("lower" = interval_lower, "upper" = interval_upper))
}

within_interval <- function(x, lower, upper){
  return(as.numeric((x >= lower) & (x <= upper)))
}

credibility_interval_selection <- function(post_sample, alpha = 0.95){
  interval <- posterior_sample_to_credibility_interval(post_sample, alpha)
  return(1 - within_interval(x = 0, lower = interval$lower, upper = interval$upper))
}







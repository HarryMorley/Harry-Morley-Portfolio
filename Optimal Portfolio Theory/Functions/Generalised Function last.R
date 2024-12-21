

###########Generalised_network_code##############

gen_network_func <- function(N, S, lambda_0, n_obs, Omega_init, tau_init, epsilons, A, catrate, mu_beta, var_beta){
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
  accrate <- array(NA, dim = c(length(epsilons), N))
  lambda_1_chain[1,,] <- exp(beta_chain[1,1] + Reduce(`+`, lapply(seq_along(A), function(j) beta_chain[j+1, 1] * A[[j]])))
  

  
  
  for(t in 1:N){
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
    
    for (i in 1:length(epsilons)){

      old_beta_likeli <- sum(dlaplace(Omega_chain[t+1, , ][lower.tri(Omega_chain[t+1, , ], diag = FALSE)], 0, exp(-(beta_chain[1, t+1] + Reduce(`+`, lapply(seq_along(A), function(j) beta_chain[j+1, t+1] * A[[j]][lower.tri(A[[j]])] )))), log = TRUE)) + dnorm(beta_chain[i, t+1], mu_beta[i], sqrt(var_beta[i]), log = TRUE)

      beta_proposal <- beta_chain[, t+1]
      beta_proposal[i] <- rnorm(1, beta_chain[i, t+1], epsilons[i])

      new_beta_likeli <- sum(dlaplace(Omega_chain[t+1, , ][lower.tri(Omega_chain[t+1, , ], diag = FALSE)], 0, exp(-(beta_proposal[1] + Reduce(`+`, lapply(seq_along(A), function(j) beta_proposal[j+1] * A[[j]][lower.tri(A[[j]])] )))), log = TRUE)) + dnorm(beta_proposal[i], mu_beta[i], sqrt(var_beta[i]), log = TRUE)

      log_alpha <- (new_beta_likeli - old_beta_likeli)

      alpha <- min(1, exp(log_alpha))

      accrate[i, t] <- alpha

      if (runif(1, 0, 1) < alpha) {
        beta_chain[i, t+1] <- beta_proposal[i]
      }
    }
    
    lambda_1_chain[t+1,,] <- exp(beta_chain[1,t+1] + Reduce(`+`, lapply(seq_along(A), function(j) beta_chain[j+1, t+1] * A[[j]])))
    
    if((t %% (N/catrate)) == 0){
      cat("Iteration", t, "done", "\n")
    }
  }
  
  
  return(list("Omega_chain" =  Omega_chain, "lambda_1_chain" =  lambda_1_chain, 'beta_chain' = beta_chain, 'accrate' = accrate))
}











############BOTH_ADJACENCY_MATRIX_CODES################




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







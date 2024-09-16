# Remote installs
#remotes::install_github("becarioprecario/INLAMSM")
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# Libraries
library(INLAMSM)
library(INLA)
library(MASS)
library(MCMCpack)
library(coda)
library(parallel)
library(BayesLogit)
library(mvnfast)
library(statmod)
library(gtools)

library(doParallel)
library(foreach)

#set.seed(123)

# Detect the number of available cores
#n_cores <- detectCores() - 1
n_cores <- 5
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Abbreviations:
# MLR : Multinomial logistic regression
# MFM : Multivariate finite mixture

# Load data set
data <- readRDS("/home/basmatsuura/Desktop/v2_thesis/R-Markdown/GibbsINLA_2/simulation/case-study-3/experiment_III-2/data.rds")
#data <- readRDS("data/synthetic/240821_univariate_spatial_n1600_noRandomY.rds")
#####################
## Global features ##

Y <- data$Y # Data set
H <- data$H # Covariates for MLR
neighbors <- dnearneigh(as.matrix(data$coords[, 1:2]), 0, 1.5)
W <- nb2mat(neighbors, style = "B")

n <- length(Y) # Number of observations
p <- 1 # Number of features
n_rho <- ncol(data$H) # Number of multinomial predictors
G <- 3 # Number of components

m <- colSums(W) # Number of neighbors per location
M <- diag(m) # Diagonal matrix for IAR precision
Q <- M-W

#########################################
##### Reference Sigma IAR effect ########

# Compute reference sigma.See Sorbey & Rue (2014)
# Sigma_star <- MASS::ginv(Q)
# diag_elements <- diag(Sigma_star)
# log_diag_elements <- log(diag_elements)
# sigma_ref <- exp(mean(0.5 * log_diag_elements))

###########################
## Gibbs-INLA parameters ##

nsim = 10
#burn = 20
#n_save <- nsim - burn

ngibbs_tot = 1000
ngibbs_burn = 800
ngibbs_save = ngibbs_tot - ngibbs_burn

#################################
####### Hyperparameters #########

gamma0 <- rep(0,n_rho) # prior mean for delta coefficients (multinomial regression)
R0k <- diag(10,n_rho) # prior covariance for delta coefficients (multinomial regression)

#########################################
######## Initialize parameters ##########

gamma_init <- matrix(0, nrow = n_rho, ncol = G-1)
phi_init <- matrix(0, nrow = n, ncol = G - 1)

eta <- cbind(rep(0,n), H%*%gamma_init + phi_init)
PI <- exp(eta) / (1 + apply(as.matrix(exp(eta[,-1])), 1, sum))

nu_init <- rep(1,G-1)

##############################
###### Helper Functions ######

psi_sums <- function(ai,Psi)
{
  Psi = as.vector(Psi)
  return(sum(Psi[ai == 1]))
}

get_psi_sums <- function(Psi, W)
{
  return(apply(W, 1, psi_sums, Psi = Psi))
}

###########################################
#### Step 1a: Update MLR parameters #######

update_mlr_params <- function(phi = phi_init, gamma = gamma_init, z = z_init, Nu2 = nu_init) {
  for (g in 1:(G-1)) {
    gamma_g <- gamma[,g]
    gamma_not_g <- gamma[,-g]
    phi_g <- phi[,g]
    phi_not_g <- phi[,-g]
    u_g <- 1*(z==(g+1))
    c_g <- log(1 + rowSums(exp(H %*% gamma_not_g + phi_not_g)))
    eta <- H %*% gamma_g + phi_g - c_g
    w <- rpg(n, 1, eta)
    u_g_star <- (u_g - 1/2)/w + c_g
    Rnk <- solve(R0k + crossprod(H*sqrt(w)))  
    rho_ng <- Rnk %*% (R0k %*% gamma0 + t(w*H) %*% u_g_star)
    gamma_g <- c(rmvn(1,rho_ng,Rnk))
    gamma_iter[,g] <- gamma_g
    
    M_g <- (get_psi_sums(phi_g, W)/m + u_g_star - H%*%gamma_g)/(m^2/Nu2[g] + 1/w^2)
    V_g <- 1/(m^2/Nu2[g] + 1/w^2)
    phi_iter[,g] <- rnorm(n = n, mean = M_g, sd = sqrt(V_g))
    Nu2[g] <- 1/rgamma(1, .01 + (n-1)/2, .01 + (t(phi_iter[,g]) %*% Q %*% phi_iter[,g])/2)
  }
  
  list(gamma = gamma_iter, phi = phi_iter, nu = Nu2)
}

#############################
####### Test Functions ######

#check_mlr_params <- update_mlr_params(phi_init, gamma_init, z_init)

###########################
### Step 1b: Fit INLA  ####

# Helper function that transforms precisions to avoid Jensen's inequality
apply_transformation_to_precision <- function(res) {
  transformed_params <- sapply(seq_len(G), function(i) {
    inla.mmarginal(
      inla.tmarginal(function(x) { 1 / (sqrt(x)) }, res$marginals.hyperpar[[i]])
    )
  })
  return(transformed_params)
}

prec.prior <- list(prec = list(prior = "logtnormal"))

fit.inla <- function(z = z_init, prec.prior = prec.prior) {
  # Initialize matrix for mixture components
  Y_mix <- matrix(NA, ncol = G, nrow = n)
  
  # Fill in the Y_mix matrix
  for (i in 1:G) {
    Y_mix[which(z == i), i] <- Y[which(z == i)]
  }
  
  # Create intercept columns
  identifier <- Y_mix
  identifier[!is.na(Y_mix)] <- 1
  
  # Manually set idx
  idx <- 1:n
  
  # Define means for fixed effects
  fixed_means <- as.list(setNames(rep(0, G), paste0("Intercept", 1:G)))
  
  # Create data list
  data_list <- list(Y_mix = Y_mix, Intercept = identifier)
  
  # Fit the INLA model
  res <<- inla(Y_mix ~ -1 + Intercept + f(idx, model = 'besag', graph = W, scale.model = TRUE),
               data = data_list,
               family = rep("gaussian", G),
               control.fixed = list(mean = fixed_means, prec = 0.01))
  
  # Transform precisions to standard deviations
  sd_vec <- apply_transformation_to_precision(res)
  
  # Reorder based on ordering of intercepts
  order_idx <- order(res$summary.fixed$mode)
  
  return(list(intercepts = res$summary.fixed$mode[order_idx],
              stand_dev = sd_vec[order_idx],
              phi_mlm = res$summary.random$idx$mode,
         nu2 = res$summary.hyperpar[G+1,]$mode))
}


check_inla <- fit.inla(data$z, prec.prior)

#####################################
### Step 2: Update latent classes ###

update_latent_classes <- function(mu, sigma, PI, phi_mlm) {
  # Step (a): Compute Pik = dnorm(yi; µk, Σk) for all observations and classes
  log_lik <- matrix(NA, nrow = n, ncol = G)
  
  for (g in 1:G) {
    log_lik[, g] <- dnorm(Y, mean = mu[g] + phi_mlm, sd = sigma[g], log = TRUE)
  }
  
  # Step (c): Compute P(zi = k|...) = Pik * πik / sum(Pih * πih) for all observations and classes
  log_Pzi_given_k <- log_lik + log(PI)
  
  # Numerical stability: subtract the maximum log value for each observation
  log_Pzi_given_k <- log_Pzi_given_k - apply(log_Pzi_given_k, 1, max)
  
  # Convert to normal scale
  Pzi_given_k <- exp(log_Pzi_given_k)
  
  # Normalize to sum to 1 for each observation
  Pzi_given_k <- Pzi_given_k / rowSums(Pzi_given_k)
  
  # Step (d): Update zi from the categorical distribution
  z_new <- apply(Pzi_given_k, 1, function(prob) sample(1:G, size = 1, prob = prob, replace = TRUE))
  
  return(z_new)
}

#######################
### Label-switching ###

remap_canonical2 <- function(z)
{
  ord_obs <- unique(z)
  z_ret <- z
  for(g in 1:length(ord_obs))
  {
    z_ret[z == ord_obs[g]] <- g
  }
  return(z_ret)
}

##################
### Gibbs-INLA ###

# Calculate the MAP estimate of the latent classes

calculate_z_map <- function(Z) {
  # Initialize a matrix to store the posterior probabilities for each class for each sample
  posterior_prob <- matrix(0, nrow = ncol(Z), ncol = G)
  
  # Calculate the posterior probability for each class
  for (i in 1:ncol(Z)) {
    posterior_prob[i, ] <- table(factor(Z[, i], levels = 1:G)) / nrow(Z)
  }
  
  # Calculate the MAP estimate
  z_map <- apply(posterior_prob, 1, which.max)
  
  return(z_map)
}

timing_result <- system.time({
  
  results <- foreach(chain = 1:n_cores, .packages = c("INLA", "BayesLogit", "mvnfast")) %dopar% {
    
    # Initialize storage for results for this chain
    chain_results <- list(
      MU = matrix(NA, nrow = nsim, ncol = G),
      SIGMA = matrix(NA, nrow = nsim, ncol = G),
      Z = matrix(NA, nrow = nsim, ncol = n),
      GAMMA = matrix(NA, nrow = nsim, ncol = (G-1) * n_rho),
      NU = matrix(NA, nrow = nsim, ncol = G-1),
      NU2 = numeric(nsim),
      Z_inner = vector("list", nsim),
      Z0 = numeric(n),
      GAMMA_inner = vector("list", nsim),
      NU_inner = vector("list", nsim)
    )
    
    # Initialize latent classes for starting Gibbs sampler in each chain
    k_means_res <- kmeans(Y, centers = G, nstart = 1)
    z_init <- k_means_res$cluster
    #z_init <- sample(1:3, size = n, replace = TRUE)
    
    # Save initial z estimate:
    chain_results$Z0 <- z_init
    
    prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01)))
    
    # Initialize variance of the IAR effects for the MLR model
    nu_init <- rep(1,G-1)
    
    # Matrices to store intermediate results in Gibbs sampler
    c_iter <- matrix(NA, nrow = n, ncol = G-1)
    eta_iter <- matrix(NA, nrow = n, ncol = G-1)
    U_star_iter <- matrix(NA, nrow = n, ncol = G-1)
    psi_sums_iter <- matrix(NA, nrow = n, ncol = G-1)
    phi_iter <- matrix(NA, nrow = n, ncol = G-1)
    M_iter <- matrix(NA, nrow = n, ncol = G-1)
    V_iter <- matrix(NA, nrow = n, ncol = G-1)
    gamma_iter <- matrix(NA, nrow = n_rho, ncol = G-1)
    
    # Gibbs sampling loop
    for (iter in 1:nsim) {
      
      inla_results <- fit.inla(z_init, prec.prior)
      mu_outer <- inla_results$intercepts
      sigma_outer <- inla_results$stand_dev
      phi_mlm_outer <- inla_results$phi_mlm
      nu2_outer <- inla_results$nu2
      
      # Storage for inner loop results
      gamma_inner <- matrix(NA, nrow = ngibbs_tot, ncol = (G-1) * n_rho)
      z_inner <- matrix(NA, nrow = ngibbs_tot, ncol = n)
      nu_inner <- matrix(NA, nrow = ngibbs_tot, ncol = G-1)
      
      #############################################
      ######## Inner Gibbs sampling loop ##########
      
      gamma_cum <- matrix(0, nrow = n_rho, ncol = G-1)
      phi_cum <- matrix(0, nrow = n, ncol = G-1)
      nu_cum <- numeric(G-1)
      z_cum <- matrix(0, nrow = ngibbs_save, ncol = n)
      
      for (g_iter in 1:ngibbs_tot) { 
        
        # Step 1: Update latent classes # MOVED UP
        z_init <- update_latent_classes(mu_outer, sigma_outer, PI, phi_mlm_outer)
        #z_init <- remap_canonical2(z_init)
        
        # Step 2: Update gamma and phi
        mlr_params <- update_mlr_params(phi_init, gamma_init, z_init, nu_init)
        gamma_init <- mlr_params$gamma
        phi_init <- mlr_params$phi
        nu_init <- mlr_params$nu
        
        # Step 3: Update multinomial logistic prior
        eta <- cbind(rep(0,n), H %*% gamma_init + phi_init)
        PI <- exp(eta) / (1 + apply(as.matrix(exp(eta[,-1])), 1, sum))
        
        # Save intermediate results from the inner loop
        gamma_inner[g_iter, ] <- c(gamma_init)
        z_inner[g_iter, ] <- z_init
        nu_inner[g_iter, ] <- nu_init
        
        if (g_iter > ngibbs_burn) {
          save_idx <- g_iter - ngibbs_burn
          
          gamma_cum <- gamma_cum + gamma_init
          phi_cum <- phi_cum + phi_init
          nu_cum <- nu_cum + nu_init
          z_cum[save_idx, ] <- z_init
        }
      }
      
      # Averaging after burn-in
      gamma_init <- gamma_cum / ngibbs_save
      phi_init <- phi_cum / ngibbs_save
      nu_init <- nu_cum / ngibbs_save
      z_init <- calculate_z_map(z_cum)
      
      # Store final results for this outer iteration
      chain_results$MU[iter, ] <- mu_outer
      chain_results$SIGMA[iter, ] <- sigma_outer
      chain_results$NU2[iter] <- nu2_outer
      chain_results$Z[iter, ] <- z_init  # Store the final z value (MAP estimate)
      chain_results$GAMMA[iter, ] <- c(gamma_init)  # Store the final gamma value
      chain_results$NU[iter, ] <- nu_init
      
      # Store the inner loop results for diagnostics, including burn-in
      chain_results$Z_inner[[iter]] <- z_inner
      chain_results$GAMMA_inner[[iter]] <- gamma_inner
      chain_results$NU_inner[[iter]] <- nu_inner
    }
    
    # Return the results for this chain
    return(chain_results)
  }
  
  # Stop the parallel cluster
  stopCluster(cl)
})

# End of Gibbs Sampler
cat("Finished Gibbs sampling after", nsim, "iterations.\n")
timing_result
save(results, data, timing_result, file = "my_data.RData")

# Aggregate results from all chains
MU_outer_list <- lapply(results, function(res) res$MU)
SIGMA_outer_list <- lapply(results, function(res) res$SIGMA)
Z_avg_list <- lapply(results, function(res) res$Z)
Z0_list <- lapply(results, function(res) res$Z0)
GAMMA_avg_list <- lapply(results, function(res) res$GAMMA)
NU_avg_list <- lapply(results, function(res) res$NU)
NU2_outer_list <- lapply(results, function(res) res$NU2)

Z_inner_list <- lapply(results, function(res) res$Z_inner)
GAMMA_inner_list <- lapply(results, function(res) res$GAMMA_inner)
NU_inner_list <- lapply(results, function(res) res$NU_inner)

##########################################################################
#### Combine Parallel Chains Results to compute Expectation of GAMMAs ####

num_matrices <- length(GAMMA_avg_list)
matrix_dim <- dim(GAMMA_avg_list[[1]])

# Initialize a matrix to store the sum of all matrices
pooled_average <- matrix(0, nrow = matrix_dim[1], ncol = matrix_dim[2])

# Sum all matrices element-wise
for (i in 1:num_matrices) {
  pooled_average <- pooled_average + GAMMA_avg_list[[i]]
}

# Divide by the number of matrices to get the element-wise average
pooled_average <- pooled_average / num_matrices

#################################################
##### Plot percenteage correctly classified #####

# Initialize data frames for Z0_list
Z0_data <- data.frame(Iteration = integer(),
                      Percentage = numeric(),
                      Chain = factor())

# Process each chain in Z0_list
for (chain in 1:length(Z0_list)) {
  if (is.vector(Z0_list[[chain]])) {
    # If it's a vector, convert it to a matrix with one row
    Z0_matrix <- matrix(Z0_list[[chain]], nrow = 1)
  } else {
    # Use it directly if it's already a matrix
    Z0_matrix <- Z0_list[[chain]]
  }
  
  # Calculate percentage of correctly classified latent classes
  percentage_matches <- mean(Z0_matrix == data$z) * 100
  
  # Append to Z0_data
  Z0_data <- rbind(Z0_data, data.frame(
    Iteration = 0,
    Percentage = percentage_matches,
    Chain = factor(chain, levels = 1:length(Z0_list))
  ))
}

# Prepare the data for Z_avg_list
plot_data <- data.frame(Iteration = rep(1:nsim, times = length(Z_avg_list)),
                        Percentage = numeric(nsim * length(Z_avg_list)),
                        Chain = factor(rep(1:length(Z_avg_list), each = nsim)))

for (chain in 1:length(Z_avg_list)) {
  percentage_matches <- apply(Z_avg_list[[chain]], 1, function(row) {
    mean(row == data$z) * 100
  })
  
  plot_data$Percentage[plot_data$Chain == chain] <- percentage_matches
}

# Combine Z0_data and plot_data
combined_plot_data <- rbind(Z0_data, plot_data)

# Create the plot using ggplot2
p <- ggplot(combined_plot_data, aes(x = Iteration, y = Percentage, color = Chain)) +
  geom_line(size = 1.2) +
  scale_y_continuous(limits = c(0, 100)) +  # Set y-axis limits from 0 to 100
  labs(title = "",
       x = "Iteration",
       y = "Percentage of Matching Latent Classes") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

# Print the plot
print(p)
ggsave("c1_matching_percentage_plot.png", plot = p, width = 8, height = 8)
##########################
###########################

# What proportion of zero-variance estimates are equivalent to the true underlying
# class?

matrix_data <- Z_inner_list[[2]][[10]]
column_variances <- apply(matrix_data, 2, var)
zero_variance_indices <- which(column_variances == 0)

zero_variance_values <- matrix_data[, zero_variance_indices]
zero_variance_values <- colMeans(zero_variance_values)

true_values <- data$z[zero_variance_indices]

matches <- zero_variance_values == true_values
proportion_matches <- mean(matches)
proportion_matches
NU_avg_list

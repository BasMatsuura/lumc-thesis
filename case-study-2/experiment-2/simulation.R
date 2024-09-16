# Libraries
library(MASS)
library(spdep)
library(mvnfast)
library(SUMMER)
library(ggplot2)

# Step 1: Define the Spatial Grid and Neighbors
dim_x <- 40
dim_y <- 40
n <- dim_x * dim_y  # Number of locations (e.g., 25x25 grid)
coords <- expand.grid(x = 1:dim_x, y = 1:dim_y)
coords$id <- 1:n

# Define neighbors using a binary proximity matrix
neighbors <- dnearneigh(as.matrix(coords[, 1:2]), 0, 1.5)
W <- nb2mat(neighbors, style = "B")

# Step 2: Generate IAR Spatial Random Effects for Multinomial Logistic Regression
G <- 5  # Number of mixture components

# Generate IAR effects (one per mixture component, except the reference component)
# phi_IAR_mlr <- matrix(NA, nrow = n, ncol = G - 1)
# for (g in 1:(G - 1)) {
#   phi_IAR_mlr[, g] <- SUMMER::rst(n = 1, type = "s", type.s = "ICAR", Amat = W, scale.model = TRUE)
#   phi_IAR_mlr <- sqrt(var_phi) * phi_IAR_mlr
# }

# Step 3: Simulate Latent Class Memberships Using the IAR Effect
p <- 4  # Number of covariates (excluding intercept)
H <- cbind(1, matrix(runif(n * p, min = -1, max = 1), n, p))  # Include intercept as the first column
gamma <- matrix(rnorm((p + 1) * (G - 1)), p + 1, G - 1)  # Coefficients for MLR (without reference category)
eta <- H %*% gamma# + phi_IAR_mlr  # Linear predictor with spatial effect
eta <- cbind(0, eta)  # Add the reference component (eta_1 = 0 for identifiability)

prob <- exp(eta) / rowSums(exp(eta))  # Multinomial probabilities
z <- apply(prob, 1, function(p) sample(1:G, 1, prob = p))  # Simulate latent classes

# Step 4: Generate IAR Spatial Random Effects for the Outcome Model

phi_IAR_mfm <- SUMMER::rst(n = 1, type = "s", type.s = "ICAR", Amat = W, scale.model = TRUE)

# Rescale unit variance IAR effect
var_phi <- 0.5  # Variance parameter for CAR model
phi_IAR_mfm <- phi_IAR_mfm * sqrt(var_phi)

# Step 4b: Generate IAR Spatial random effect using joint distribution
# diag(colSums(W))[1,1]
# Q <- diag(colSums(W)) - W
# cov_mat <- solve(Q + diag(0.000001, n)) 
# phi_IAR_mfm_joint <- mvnfast::rmvn(1, rep(0,n), cov_mat)

# Step 5: Generate Multivariate Observations Y with mCAR Effect

# Hardcoded means and standard deviations
mu <- c(-7, -3.5, 0, 3.5, 7)
Sigma <- c(1,1,1,1,1)

Y <- numeric(n)
for (i in 1:n) {
  g <- z[i]
  Y[i] <- rnorm(1, mu[g] + phi_IAR_mfm[i], Sigma[g])  # Add mCAR effect to the mean
}

# The synthetic dataset
synthetic_data <- list(
  coords = coords,
  W = W,
  z = z,
  Y = Y,
  H = H,
  #phi_IAR_mlr = phi_IAR_mlr,
  phi_IAR_mfm = phi_IAR_mfm,
  mu = mu,
  Sigma = Sigma,
  gamma = gamma
)

# Reshape one of the phi_IAR effects and phi_mCAR to be matrices for the image function
#phi_IAR_mlr_matrix <- matrix(phi_IAR_mlr[, 1], nrow = dim_x, ncol = dim_y)
#phi_IAR_mfm_matrix <- matrix(phi_IAR_mfm[, 1], nrow = dim_x, ncol = dim_y)

# Plotting the results
#par(mfrow = c(1, 3))
#plot(coords$x, coords$y, col = z, pch = 19, main = "Latent Classes (z)")
#image(1:dim_x, 1:dim_y, t(phi_IAR_mlr_matrix)[, ncol(phi_IAR_mlr_matrix):1], main = "IAR Effect MLR (Component 1)", col = heat.colors(n))
#image(1:dim_x, 1:dim_y, t(phi_IAR_mfm_matrix)[, ncol(phi_IAR_mfm_matrix):1], main = "IAR Effect MFM", col = heat.colors(n))

# Define the number of bins
num_bins <- 50  # Adjust this value as needed

# Create a data frame for plotting
plot_data <- data.frame(Y = synthetic_data$Y, Cluster = as.factor(synthetic_data$z))

# Create the histogram plot with the specified number of bins
p_histogram <- ggplot(plot_data, aes(x = Y, fill = Cluster)) +
  geom_histogram(bins = num_bins, position = "identity", alpha = 0.6) +  # Use the num_bins parameter
  scale_fill_viridis_d(name = "Cluster") +  # Use a color palette for clusters
  labs(title = "",
       x = "Data Value",
       y = "Frequency") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Print the plot
print(p_histogram)

table(synthetic_data$z)
# Save file
#saveRDS(synthetic_data, file = "/home/basmatsuura/Desktop/v2_thesis/R-Markdown/GibbsINLA_2/synthetic_data.rds")


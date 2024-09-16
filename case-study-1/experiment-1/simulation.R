# Libraries
library(MASS)
library(spdep)
library(mvnfast)
#library(SUMMER)

# Set seed
#set.seed(123)

n <- 1500

G <- 10  # Number of mixture components

p <- 7 
H <- cbind(1, matrix(rnorm(n * p), n, p)) 
gamma <- matrix(rnorm((p + 1) * (G - 1)), p + 1, G - 1)  #
eta <- H %*% gamma#
eta <- cbind(0, eta)  
prob <- exp(eta) / rowSums(exp(eta))  # Multinomial probabilities
z <- apply(prob, 1, function(p) sample(1:G, 1, prob = p))  # Simulate latent classes


mu <- c(-9,-7,-5,-3,-1,1,3,5,7,9)
Sigma <- rep(1,10)

Y <- numeric(n)
for (i in 1:n) {
  g <- z[i]
  Y[i] <- rnorm(1, mu[g] , Sigma[g])  # Add mCAR effect to the mean
}

# The synthetic dataset
synthetic_data <- list(
  coords = coords,
  z = z,
  Y = Y,
  H = H,
  mu = mu,
  Sigma = Sigma,
  gamma = gamma
)

# Define the number of bins
num_bins <- 30  # Adjust this value as needed

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
saveRDS(synthetic_data, file = "/home/basmatsuura/Desktop/v2_thesis/R-Markdown/GibbsINLA_2/synthetic_data.rds")


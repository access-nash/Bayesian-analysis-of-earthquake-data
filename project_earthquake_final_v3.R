# Load required libraries
library(MASS)        # For mvrnorm
library(MCMCpack)    # For rdirichlet
library(ggplot2)     # For plotting

# Load earthquake data
earthquakes.dat <- read.delim("C:/Users/avina/OneDrive/Documents/earthquakes.txt")
earthquakes.dat$Quakes <- as.numeric(earthquakes.dat$Quakes)
y.dat <- earthquakes.dat$Quakes[1:100]  # Training data

# Exploratory plot of the data
ggplot(earthquakes.dat, aes(x = 1:103, y = Quakes)) +
  geom_line() +
  labs(title = "Earthquake Data Time Series", x = "Time", y = "Quakes") +
  theme_minimal()

# Autocorrelation and Partial Autocorrelation Plots to check for stationarity and lags
acf(y.dat, main = "ACF of Earthquake Data")
pacf(y.dat, main = "PACF of Earthquake Data")

# Define necessary parameters
num_iterations <- 10000
burn_in <- 5000
num_components <- 2  # For mixture model
p <- 3  # AR(3) model order

# Prepare lagged matrix X for AR(3)
n <- length(y.dat)
X <- cbind(y.dat[3:(n-1)], y.dat[2:(n-2)], y.dat[1:(n-3)])
y <- y.dat[4:n]

# Set weakly informative priors
a1 <- 1
a2 <- 1
m0 <- rep(0, p)
C0 <- diag(10, p)  # This is the missing covariance matrix
n0 <- 0.02
d0 <- 0.02

# Gibbs sampler for AR(3) model
gibbs_sampler_ar <- function(y, X, num_iterations, burn_in, m0, C0) {
  n <- nrow(X)
  p <- ncol(X)
  phi_samples <- matrix(0, nrow = num_iterations, ncol = p)
  sigma2_samples <- numeric(num_iterations)
  
  # Initial values
  phi <- rep(0, p)
  sigma2 <- 1
  
  for (iter in 1:num_iterations) {
    # Update phi
    Sigma_phi <- solve(t(X) %*% X / sigma2 + solve(C0))
    mu_phi <- Sigma_phi %*% (t(X) %*% y / sigma2 + solve(C0) %*% m0)
    phi <- mvrnorm(1, mu = mu_phi, Sigma = Sigma_phi)
    
    # Update sigma2
    residuals <- y - X %*% phi
    shape <- n0 + n / 2
    scale <- d0 + sum(residuals^2) / 2
    sigma2 <- 1 / rgamma(1, shape, rate = scale)
    
    # Store samples
    phi_samples[iter, ] <- phi
    sigma2_samples[iter] <- sigma2
  }
  
  # Burn-in
  phi_samples <- phi_samples[(burn_in+1):num_iterations, ]
  sigma2_samples <- sigma2_samples[(burn_in+1):num_iterations]
  
  return(list(phi = phi_samples, sigma2 = sigma2_samples))
}

# Fit AR(3) model
posterior_samples_ar3 <- gibbs_sampler_ar(y, X, num_iterations, burn_in, m0, C0)

# Posterior analysis: Plot posterior distributions of the coefficients
par(mfrow = c(2, 2))
for (i in 1:3) {
  hist(posterior_samples_ar3$phi[, i], main = paste("Posterior of phi", i), xlab = paste("phi", i), probability = TRUE)
}
hist(posterior_samples_ar3$sigma2, main = "Posterior of sigma^2", xlab = "sigma^2", probability = TRUE)


# Gibbs sampler for mixture of AR(3) model with checks for empty components
gibbs_sampler_mixture_ar <- function(y, X, num_components, num_iterations, burn_in, m0, C0) {
  n <- nrow(X)
  p <- ncol(X)
  phi_samples <- array(0, dim = c(num_iterations, num_components, p))
  sigma2_samples <- matrix(0, nrow = num_iterations, ncol = num_components)
  pi_samples <- matrix(0, nrow = num_iterations, ncol = num_components)
  
  # Initialize values
  phi <- matrix(0, nrow = num_components, ncol = p)
  sigma2 <- rep(1, num_components)
  pi_k <- rep(1 / num_components, num_components)
  z <- sample(1:num_components, n, replace = TRUE)
  
  alpha_dirichlet <- rep(1, num_components)
  
  for (iter in 1:num_iterations) {
    # Update component assignments z
    for (i in 1:n) {
      probs <- pi_k * dnorm(y[i], mean = X[i, ] %*% t(phi), sd = sqrt(sigma2))
      z[i] <- sample(1:num_components, 1, prob = probs / sum(probs))
    }
    
    # Update component probabilities pi_k
    nk <- table(factor(z, levels = 1:num_components))
    pi_k <- rdirichlet(1, alpha_dirichlet + nk)
    
    # Update phi and sigma2 for each component
    for (k in 1:num_components) {
      X_k <- X[z == k, , drop = FALSE]  # Select the rows of X assigned to component k
      y_k <- y[z == k]
      
      if (nrow(X_k) > 0) {  # Check if X_k is not empty
        Sigma_phi <- solve(t(X_k) %*% X_k / sigma2[k] + solve(C0))
        mu_phi <- Sigma_phi %*% (t(X_k) %*% y_k / sigma2[k] + solve(C0) %*% m0)
        phi[k, ] <- mvrnorm(1, mu = mu_phi, Sigma = Sigma_phi)
        
        residuals <- y_k - X_k %*% phi[k, ]
        shape <- n0 + length(y_k) / 2
        scale <- d0 + sum(residuals^2) / 2
        sigma2[k] <- 1 / rgamma(1, shape, rate = scale)
      }
    }
    
    # Store samples
    phi_samples[iter, , ] <- phi
    sigma2_samples[iter, ] <- sigma2
    pi_samples[iter, ] <- pi_k
  }
  
  # Burn-in
  phi_samples <- phi_samples[(burn_in+1):num_iterations, , ]
  sigma2_samples <- sigma2_samples[(burn_in+1):num_iterations, ]
  pi_samples <- pi_samples[(burn_in+1):num_iterations, ]
  
  return(list(phi = phi_samples, sigma2 = sigma2_samples, pi = pi_samples))
}

# Fit mixture AR(3) model
posterior_samples_mixture <- gibbs_sampler_mixture_ar(y, X, num_components, num_iterations, burn_in, m0, C0)

# Posterior summary of AR(3)
summary_phi_ar3 <- apply(posterior_samples_ar3$phi, 2, mean)
summary_sigma2_ar3 <- mean(posterior_samples_ar3$sigma2)
cat("Posterior Mean of AR(3) Coefficients:", summary_phi_ar3, "\n")
cat("Posterior Mean of AR(3) Variance:", summary_sigma2_ar3, "\n")

# Predictions for AR(3) model (3-step ahead)
X_pred <- cbind(y.dat[98:100], y.dat[97:99], y.dat[96:98])
predictions_ar3 <- X_pred %*% summary_phi_ar3
cat("3-step ahead predictions from AR(3):", predictions_ar3, "\n")

# Function to compute AIC, BIC, and DIC for AR models
compute_aic_bic_dic <- function(y, X, posterior_samples, p, num_components) {
  # Calculate the number of parameters
  n_params <- p * num_components + num_components  # AR coefficients and sigma2 for each component
  
  # Extract the posterior mean of phi and sigma2
  # Assuming posterior_samples$phi is a matrix of dimension (n_iterations x p)
  phi_mean <- apply(posterior_samples$phi, 2, mean)  # Mean of phi across iterations
  sigma2_mean <- mean(posterior_samples$sigma2)  # Mean of sigma2 across iterations
  
  # Check dimensions of phi_mean and X
  if (ncol(X) != length(phi_mean)) {
    stop("Dimensions of X and phi_mean do not match.")
  }
  
  # Compute log-likelihood
  log_likelihood <- 0
  for (i in 1:nrow(X)) {
    mu_i <- X[i, ] %*% phi_mean  # Linear prediction
    log_likelihood <- log_likelihood + dnorm(y[i], mean = mu_i, sd = sqrt(sigma2_mean), log = TRUE)
  }
  
  # AIC: -2 * log-likelihood + 2 * number of parameters
  aic <- -2 * log_likelihood + 2 * n_params
  
  # BIC: -2 * log-likelihood + log(n) * number of parameters
  bic <- -2 * log_likelihood + log(length(y)) * n_params
  
  # DIC: Mean of the deviance minus deviance at posterior mean
  deviance_post_mean <- -2 * sum(dnorm(y, mean = X %*% phi_mean, sd = sqrt(sigma2_mean), log = TRUE))
  deviance_samples <- apply(posterior_samples$phi, 1, function(phi_sample) {
    -2 * sum(dnorm(y, mean = X %*% phi_sample, sd = sqrt(posterior_samples$sigma2), log = TRUE))
  })
  dic <- mean(deviance_samples) + (mean(deviance_samples) - deviance_post_mean)
  
  return(list(aic = aic, bic = bic, dic = dic))
}

# Use the function to compute AIC, BIC, DIC for AR(3) model
ar3_criteria <- compute_aic_bic_dic(y, X, posterior_samples_ar3, p = 3, num_components = 1)

# Print the criteria
print(ar3_criteria)

# Function to compute AIC, BIC, and DIC for Mixture AR models
compute_aic_bic_dic_mixture <- function(y, X, posterior_samples, p, num_components) {
  # Calculate the number of parameters
  n_params <- p * num_components + num_components  # AR coefficients and sigma2 for each component
  
  # Extract posterior means for each component
  phi_means <- list()
  sigma2_means <- numeric(num_components)
  for (k in 1:num_components) {
    phi_means[[k]] <- apply(posterior_samples$phi[, k, ], 2, mean)  # Mean of phi for each component
    sigma2_means[k] <- mean(posterior_samples$sigma2[, k])  # Mean of sigma2 for each component
  }
  
  # Check dimensions of X and phi_means for each component
  for (k in 1:num_components) {
    if (ncol(X) != length(phi_means[[k]])) {
      stop("Dimensions of X and phi_mean for component ", k, " do not match.")
    }
  }
  
  # Calculate the log-likelihood for the mixture model
  log_likelihood <- 0
  for (i in 1:nrow(X)) {
    log_likelihood_i <- 0
    for (k in 1:num_components) {
      mu_i <- X[i, ] %*% phi_means[[k]]  # Linear prediction for component k
      log_likelihood_i <- log_likelihood_i + dnorm(y[i], mean = mu_i, sd = sqrt(sigma2_means[k]), log = TRUE)
    }
    log_likelihood <- log_likelihood + log_likelihood_i
  }
  
  # AIC: -2 * log-likelihood + 2 * number of parameters
  aic <- -2 * log_likelihood + 2 * n_params
  
  # BIC: -2 * log-likelihood + log(n) * number of parameters
  bic <- -2 * log_likelihood + log(length(y)) * n_params
  
  # DIC: Mean of the deviance minus deviance at posterior mean
  deviance_post_mean <- 0
  for (k in 1:num_components) {
    deviance_post_mean <- deviance_post_mean -2 * sum(dnorm(y, mean = X %*% phi_means[[k]], sd = sqrt(sigma2_means[k]), log = TRUE))
  }
  
  deviance_samples <- numeric(nrow(posterior_samples$phi))
  for (i in 1:nrow(posterior_samples$phi)) {
    deviance_i <- 0
    for (k in 1:num_components) {
      phi_sample <- posterior_samples$phi[i, k, ]
      deviance_i <- deviance_i -2 * sum(dnorm(y, mean = X %*% phi_sample, sd = sqrt(posterior_samples$sigma2[i, k]), log = TRUE))
    }
    deviance_samples[i] <- deviance_i
  }
  dic <- mean(deviance_samples) + (mean(deviance_samples) - deviance_post_mean)
  
  return(list(aic = aic, bic = bic, dic = dic))
}

# Use the function to compute AIC, BIC, DIC for the Mixture AR(3) model
mixture_criteria <- compute_aic_bic_dic_mixture(y, X, posterior_samples_mixture, p = 3, num_components = 2)

# Print the criteria for the Mixture AR model
print(mixture_criteria)

# Posterior summary of Mixture AR(3)
summary_phi_mixture <- apply(posterior_samples_mixture$phi, c(2, 3), mean)
summary_sigma2_mixture <- apply(posterior_samples_mixture$sigma2, 2, mean)
summary_lambda_mixture <- mean(posterior_samples_mixture$lambda)
cat("Posterior Mean of Mixture AR(3) Coefficients:", summary_phi_mixture, "\n")
cat("Posterior Mean of Mixture AR(3) Variance:", summary_sigma2_mixture, "\n")
cat("Posterior Mean of Lambda (Mixing proportion):", summary_lambda_mixture, "\n")

# Predictions for Mixture AR(3) model (3-step ahead)
predictions_mixture <- X_pred %*% summary_phi_mixture[, 1] * summary_lambda_mixture +
  X_pred %*% summary_phi_mixture[, 2] * (1 - summary_lambda_mixture)
cat("3-step ahead predictions from Mixture AR(3):", predictions_mixture, "\n")


# Compare models
print(ar3_criteria)
print(mixture_criteria)

# Conclusion
cat("AR(3) Model AIC, BIC, DIC:", ar3_criteria, "\n")
cat("Mixture AR(3) Model AIC, BIC, DIC:", mixture_criteria, "\n")

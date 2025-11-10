
# NOTE - this method looks nice, but is not strictly the process we want. By
# modelling a gaussian field in 3D we find that probability can "bleed" through
# the sphere, so that points on opposite sides are only 2*r apart, and not pi*r
# as we would expect

# Load necessary libraries
library(MASS) # For mvrnorm function to generate multivariate normal samples
library(pracma) # For dot product and norm calculations
library(glmnet)
library(plotly)

# Function to sample points uniformly on a sphere
sample_sphere <- function(n, R) {
  phi <- runif(n, 0, 2 * pi)
  costheta <- runif(n, -1, 1)
  theta <- acos(costheta)
  
  x <- R * sin(theta) * cos(phi)
  y <- R * sin(theta) * sin(phi)
  z <- R * cos(theta)
  
  return(cbind(x, y, z))
}

# Function to generate random Fourier features with lengthscale parameter
generate_rff <- function(n_features, R, lengthscale) {
  w <- sample_sphere(n_features, R) / lengthscale
  return(list(w = w))
}

# Function to compute the feature map
feature_map <- function(x, rff) {
  n_features <- nrow(rff$w)
  z <- matrix(0, nrow = n_features, ncol = 2)
  for (i in 1:n_features) {
    w_dot_x <- sum(rff$w[i, ] * x)
    z[i, 1] <- cos(w_dot_x)
    z[i, 2] <- sin(w_dot_x)
  }
  z <- sqrt(2 / n_features) * as.vector(t(z))
  return(z)
}

# ------------------------------------------------

# Set parameters
R <- 10 # Radius of the sphere
n_features <- 100 # Number of random features
n_points <- 500 # Number of points on the sphere
lengthscale <- 20 # note this interacts with the radius parameter

# Generate random points on the sphere (for demonstration purposes)
#set.seed(123)
points <- sample_sphere(n_points, R)

# Generate random Fourier features
rff <- generate_rff(n_features, R, lengthscale)

# Compute feature maps for all points
feature_maps <- t(apply(points, 1, function(p) feature_map(p, rff)))

# Generate some synthetic data for demonstration purposes
true_alpha <- rnorm(n_features * 2)
y <- apply(feature_maps, 1, function(z) sum(true_alpha * z))# + rnorm(n_points, 0, 0.1)

plot_ly(x = ~x, y = ~y, z = ~z, color = ~q, type = 'scatter3d', mode = 'markers',
        marker = list(size = 5),
        data = data.frame(x = points[,"x"],
                          y = points[,"y"],
                          z = points[,"z"],
                          q = y))

# ------------------------------------------------

# Fit Ridge regression model
fit_ridge <- glmnet(feature_maps, y, alpha = 0, intercept = TRUE) # alpha = 0 for Ridge regression

# Get the coefficients
estimated_alpha_ridge <- coef(fit_ridge, s = 0.01)[-1,1] # 's' is the regularization parameter

# ------------------------------------------------

# Function to predict values using the fitted model
predict_grf <- function(new_points, rff, estimated_alpha) {
  new_feature_maps <- t(apply(new_points, 1, function(p) feature_map(p, rff)))
  predicted_values <- new_feature_maps %*% estimated_alpha
  return(predicted_values)
}

# Generate some new points on the sphere
new_points <- sample_sphere(1e4, R)

# Predict values at new points
predicted_values <- predict_grf(new_points, rff, estimated_alpha_ridge)

plot_ly(x = ~x, y = ~y, z = ~z, color = ~q, type = 'scatter3d', mode = 'markers',
        marker = list(size = 5),
        data = data.frame(x = new_points[,"x"],
                          y = new_points[,"y"],
                          z = new_points[,"z"],
                          q = predicted_values))


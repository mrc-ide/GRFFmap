# 11.length_scale_check.R
#
# Author: Bob Verity
# Date: 2024-08-02
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# Checks that all three simulation methods produce the same length-scale.
# Methods: 1) spectral, 2) naive mvnorm draw, 3) RFF
#
# Conclusion: don't forget to include the factor of 1/2 in the RBF kernel!
#
# ------------------------------------------------------------------

library(mvtnorm)

# parameters
x_range <- c(-180, 180)
y_range <- c(-90, 90)
t_range <- c(1, 1)
nx <- 50
ny <- 40
length_space <- 5
sill <- 3
kernel_space_type <- "RBF" # options {RBF, Exp}
kernel_time_type <- "RBF"


# simulate grid
sim_array <- sim_spacetime_spectral(x_range = x_range,
                                    y_range = y_range,
                                    t_range = t_range,
                                    nx = nx,
                                    ny = ny,
                                    length_space = length_space,
                                    length_time = 1,
                                    x_buffer = nx / 2,
                                    y_buffer = ny / 2,
                                    t_buffer = 0,
                                    kernel_space_type = kernel_space_type,
                                    kernel_time_type = kernel_time_type,
                                    sill = sill)

m1 <- sim_array$grid[,,1]
image(m1)

# ------------------------------------------------

df_coords <- expand_grid(x = sim_array$x,
                         y = sim_array$y)

d <- dist(df_coords) |>
  as.matrix()
rho <- exp(-0.5*(d / length_space)^2)

z_draw <- rmvnorm(1, sigma = sill*rho)[1,]
m2 <- t(matrix(z_draw, nrow = ny))
image(m2)

# ------------------------------------------------

n_features <- 100
omega_space <- rbind(rnorm(n_features),
                     rnorm(n_features))
X <- as.matrix(df_coords)
feat <- cbind(cos(X %*% omega_space / length_space),
              sin(X %*% omega_space / length_space)) / sqrt(n_features)

beta <- rnorm(n_features * 2, sd = sqrt(sill))
z_draw2 <- as.vector(feat %*% beta)

m3 <- t(matrix(z_draw2, nrow = ny))
image(m3)

# ------------------------------------------------

par(mfrow = c(2,2))
image(m1, zlim = c(-5, 5))
image(m2, zlim = c(-5, 5))
image(m3, zlim = c(-5, 5))

var(as.vector(m1))
var(as.vector(m2))
var(as.vector(m3))



library(tidyverse)
library(plotly)
library(rstan)
library(av)

get_fade <- function(t_vec, t_data, show_window, fade_window) {
  
  # calculate fade matrix
  d <- outer(t_vec, t_data, FUN = "-")
  ret <- 1 - (d - show_window) / fade_window
  ret[ret > 1] <- 1
  ret[ret < 0] <- 0
  ret[d < 0] <- 0
  
  ret
}


# parameters
x_range <- c(-180, 180)
y_range <- c(-90, 90)
t_range <- c(1, 100)
nx <- 120
ny <- 80
length_space <- 30
length_time <- 20
sill <- 2
kernel_space_type <- "RBF" # options {RBF, Exp}
kernel_time_type <- "RBF"


# simulate grid
sim_array <- sim_spacetime_spectral(x_range = x_range,
                                    y_range = y_range,
                                    t_range = t_range,
                                    nx = nx,
                                    ny = ny,
                                    length_space = length_space,
                                    length_time = length_time,
                                    x_buffer = 50,
                                    y_buffer = 50,
                                    t_buffer = 50,
                                    kernel_space_type = kernel_space_type,
                                    kernel_time_type = kernel_time_type,
                                    sill = sill)

image(sim_array$grid[,,1])

image(sim_array$grid[1,,])

hist(as.vector(sim_array$grid))
sd(as.vector(sim_array$grid))

# simulate data from grid
n_data <- 400
df_dat <- array_draw_data(sim_array = sim_array,
                          n_data = n_data,
                          N = 10)

# get distance between observations in space, time, and value
data_dist <- get_data_dist(x = df_dat$x,
                           y = df_dat$y,
                           t = df_dat$t,
                           z = df_dat$p_emplogit)

# plot_ly(data_dist, x = ~dist_space, y = ~dist_time, z = ~dist_z, size = 1,
#         type = "scatter3d", mode = "markers")

if (FALSE) {
  # get loglikelihood over all parameter combinations
  df_param_combos <- expand_grid(length_space = seq(1, 100, l = 21),
                                 length_time = seq(1, 100, l = 21),
                                 partial_sill = seq(0, 3, l = 21),
                                 nugget = seq(0, 1, l = 21))
  df_param_combos$ll <- mapply(function(i) {
    if ((i %% 1e3) == 0) {
      message(sprintf("%s of %s", i, nrow(df_param_combos)))
    }
    get_vario_loglike(dist_space = data_dist$dist_space,
                      dist_time = data_dist$dist_time,
                      dist_z = data_dist$dist_z,
                      length_space = df_param_combos$length_space[i],
                      length_time = df_param_combos$length_time[i],
                      nugget = df_param_combos$nugget[i],
                      partial_sill = df_param_combos$partial_sill[i],
                      kernel_space_power = 1,
                      kernel_time_power = 2)
  }, 1:nrow(df_param_combos))
  
  # extract ML values
  w <- which.max(df_param_combos$ll)
  length_space_est <- df_param_combos$length_space[w]
  length_time_est <- df_param_combos$length_time[w]
  nugget_est <- df_param_combos$nugget[w]
  partial_sill_est <- df_param_combos$partial_sill[w]
  
  # plot likelihood surface
  df_param_combos |>
    filter(length_space == length_space_est) |>
    filter(length_time == length_time_est) |>
    mutate(like = exp(ll - max(ll, na.rm = TRUE))) |>
    ggplot() + theme_bw() +
    geom_raster(aes(x = partial_sill, y = nugget, fill = like))
  
  df_param_combos |>
    filter(partial_sill == partial_sill_est) |>
    filter(nugget == nugget_est) |>
    mutate(like = exp(ll - max(ll, na.rm = TRUE))) |>
    ggplot() + theme_bw() +
    geom_raster(aes(x = length_space, y = length_time, fill = like))
  
  length_space
  length_space_est
  
  length_time
  length_time_est
}

# ------------------------------------------------

# draw features
n_features <- 50
if (kernel_space_type == "RBF") {
  omega_space <- rbind(rnorm(n_features),
                       rnorm(n_features))
} else if (kernel_space_type == "Exp") {
  omega_space <- rbind(rcauchy(n_features),
                       rcauchy(n_features))
}
omega_time <- rbind(rnorm(n_features))

# make feature maps
X <- df_dat |>
  select(x, y) |>
  as.matrix()

feat <- cbind(cos(X %*% omega_space / length_space + df_dat$t %*% omega_time / length_time),
              sin(X %*% omega_space / length_space + df_dat$t %*% omega_time / length_time)) / sqrt(n_features)

# Prepare data list for Stan
data_list <- list(
  N = n_data,
  n_features = n_features,
  feat = feat,
  k = df_dat$k,
  n_trials = df_dat$N,
  psill_shape = 1,
  psill_rate = 1 / 100,
  nugget_shape = 1,
  nugget_rate = 1,
  mu_mean = 0,
  mu_sd = 1
)

# run stan model
fit <- stan(
  file = "R_ignore/R_scripts/stan/fixedlength_vario.stan",
  data = data_list,
  iter = 1e3,
  chains = 1
)

# Extract samples for beta and std
draws <- extract(fit)
beta_samples <- draws$beta
sill_samples <- draws$psill
n_draws <- nrow(beta_samples)

hist(sill_samples)

# ------------------------------------------------

# make feature maps for prediction
X_pred <- expand.grid(x = sim_array$x,
                      y = sim_array$y) |>
  as.matrix()

# get fade
t_vec <- sim_array$t
fade_mat <- get_fade(t_vec, df_dat$t, show_window = 10, fade_window = 10)

for (i in seq_along(t_vec)) {
  message(sprintf("%s of %s", i, length(t_vec)))
  
  # fix prediction at this time
  t_pred <- rep(t_vec[i], nrow(X_pred))
  
  # generate prediction feature maps
  feat_pred <- cbind(cos(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time),
                     sin(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time)) / sqrt(n_features)
  
  # get predictions over all posterior draws
  post_all <- feat_pred %*% t(beta_samples) |>
    plogis()
  
  # summarise posterior predictions
  pred_quant <- apply(post_all, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) |>
    t() |>
    as.data.frame()
  
  # get plotting data.frames from infered distribution and truth
  df_pred <- data.frame(x = X_pred[,1],
                        y = X_pred[,2],
                        t = t_pred) |>
    bind_cols(pred_quant) |>
    mutate(MOE = `97.5%` - `2.5%`,
           p_pred = `50%`,
           mask = (MOE > 10.25))
  
  df_true <- data.frame(x = X_pred[,1],
                        y = X_pred[,2],
                        t = t_pred,
                        z = as.vector(sim_array$grid[,,t_vec[i]])) |>
    mutate(p_true = plogis(z))
  
  # produce plots
  plot1 <- df_true |>
    ggplot() + theme_bw() +
    geom_raster(aes(x = x, y = y, fill = p_true)) +
    geom_point(aes(x = x, y = y, fill = p_est, alpha = alpha),
               size = 3, pch = 21, data = cbind(df_dat, alpha = fade_mat[i,])) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
    theme(legend.position = "none") +
    ggtitle("True prevalence")
  #plot1
  
  df_masked <- df_pred |>
    filter(mask == TRUE)
  dx <- diff(sim_array$x)[1]
  dy <- diff(sim_array$y)[1]
  plot2 <- df_pred |>
    ggplot() + theme_bw() +
    geom_raster(aes(x = x, y = y, fill = p_pred)) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
    geom_tile(aes(x = x, y = y), width = dx, height = dy, fill = grey(0.5), alpha = 0.4, data = df_masked) +
    theme(legend.position = "none") +
    ggtitle("Inferred prevalence")
  #plot2
  
  # combine plots
  plot_combined <- cowplot::plot_grid(plot1, plot2)
  #plot_combined
  
  # save to file as .png
  ggsave(sprintf("R_ignore/R_scripts/outputs/sim2_frames/frame%s.png", i),
         plot = plot_combined, width = 8, height = 4, dpi = 150)
  
}

image_files <- sprintf("R_ignore/R_scripts/outputs/sim2_frames/frame%s.png", seq_along(t_vec))
av::av_encode_video(input = image_files,
                    output = "R_ignore/R_scripts/outputs/sim2.mp4",
                    framerate = 20)

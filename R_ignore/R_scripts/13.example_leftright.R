# 13.example_leftright.R
#
# Author: Bob Verity
# Date: 2024-09-25
#
# Inputs: (none)
#
# Outputs: /Users/rverity/Desktop/GRFF_example1.rds
#
# Purpose:
# Tests that the model gets the correct inference given a known prevalence
# surface changing over time. Create a known prevalence space-time grid in which
# prevalence is initially high in the west but moves to the east over time.
# Saves the output in a list by time slice allowing downstream mapping.
#
# ------------------------------------------------------------------

# define model parameters
n_sites <- 100
N <- 100

lat_min <- -0.87
lat_max <- 3.30
lon_min <- 31.38
lon_max <- 33.77
t_min <- 1
t_max <- 7000
nx <- 50
ny <- nx
nt <- 51

# make a cube of true prevalence
df_p_true <- expand_grid(lat = seq(0, 1, l = ny),
                         lon = seq(0, 1, l = nx),
                         t = seq(0, 1, l = nt)) |>
  mutate(p_true = (1 - lon)^2*(1 - t) + lon^2*t) |>
  #mutate(p_true = (1 - lon)^2) |>
  mutate(lat = lat_min + lat*(lat_max - lat_min),
         lon = lon_min + lon*(lon_max - lon_min),
         t = round(t_min + t*(t_max - t_min)))

# plot true prevalence at time slice
t_vec <- sort(unique(df_p_true$t))
t_plot <- t_vec[1]
df_p_true |>
  filter(t == t_plot) |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = lon, y = lat, fill = p_true)) +
  scale_fill_viridis_c(option = "magma")

# draw data at random sites
prev_data <- df_p_true |>
  sample_n(n_sites) |>
  mutate(denominator = N) |>
  mutate(numerator = rbinom(n = n(), size = denominator, prob = p_true), .before = denominator)

head(prev_data)

# ------------------------------------------------

# prepare for model fitting
prev_data <- prev_data |>
  mutate(p_emp = emplogit(numerator, denominator))

# plot data
prev_data |>
  ggplot() + theme_bw() +
  geom_point(aes(x = lon, y = lat, size = t, col = p_emp))

# define model parameters
length_space <- 1
length_time <- 365
n_features <- 50

# draw features
omega_space <- rbind(rnorm(n_features),
                     rnorm(n_features))
omega_time <- rbind(rnorm(n_features))

# make feature maps
X <- prev_data |>
  dplyr::select(x = lon, y = lat) |>
  as.matrix()

feat <- cbind(cos(X %*% omega_space / length_space + prev_data$t %*% omega_time / length_time),
              sin(X %*% omega_space / length_space + prev_data$t %*% omega_time / length_time)) / sqrt(n_features)

# Prepare data list for Stan
data_list <- list(
  n_data = n_sites,
  n_features = n_features,
  feat = feat,
  n = prev_data$denominator,
  z = prev_data$p_emp
)

# ------------------------------------------------

# run stan model
fit <- stan(
  file = "R_ignore/R_scripts/stan/logistic_weighted.stan",
  data = data_list,
  iter = 1e4,
  chains = 1
)

summary_stats <- summary(fit)$summary
summary_stats[,c("n_eff", "Rhat")]

# extract draws
thinning <- 10

draws <- rstan::extract(fit)
thin_vec <- seq(1, length(draws$sigma), by = thinning)
n_draws_thin <- length(thin_vec)

beta_samples <- draws$beta[thin_vec,]
sigma_samples <- draws$sigma[thin_vec]
mu_samples <- draws$mu[thin_vec]

# ------------------------------------------------

# read in prediction grid
r <- raster("/Users/rverity/Dropbox/Bob/Work/Events/2024.09.23 IDEEL hackathon London/analysis/stitch/analysis/data_derived/uganda_raster_admin0.tif")
X_pred <- coordinates(r)
n_pred <- nrow(X_pred)

# set times at which to predict
t_vec <- seq(1, 7000, 30)

# loop over time
l <- list()
for (i in seq_along(t_vec)) {
  message(sprintf("time %s of %s", i, length(t_vec)))
  
  # create fixed time vector
  t_pred <- rep(t_vec[i], n_pred)
  
  # generate prediction feature maps
  feat_pred <- cbind(cos(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time),
                     sin(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time)) / sqrt(n_features)
  
  # get predictions over all posterior draws
  post_all <- feat_pred %*% t(beta_samples) |>
    sweep(2, mu_samples, "+") |>
    plogis()
  
  # get a series of posterior draws to output
  post_draws <- post_all[,sample(n_draws_thin, 100)] |>
    as.data.frame() |>
    mutate(lat = X_pred[,2],
           lon = X_pred[,1], .before = 1)
  
  # summarise posterior predictions
  post_summary <- apply(post_all, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) |>
    t() |>
    as.data.frame() |>
    rename(prev_Q2.5 = '2.5%',
           prev_Q50 = '50%',
           prev_Q97.5 = '97.5%') |>
    mutate(prev_sd = apply(X = post_all, MARGIN = 1, FUN = sd)) |>
    mutate(lat = X_pred[,2],
           lon = X_pred[,1], .before = 1)
  
  post_summary |>
    ggplot() + theme_void() +
    geom_raster(aes(x = lon, y = lat, fill = prev_Q50)) +
    scale_fill_viridis_c(option = "plasma", limits = c(0, 1))
  
  l[[i]] <- list(time = t_vec[i],
               post_summary = post_summary,
               post_draws = post_draws)
}

bobfunctions2::object.size_auto(l)

if (FALSE) {
  saveRDS(l, "/Users/rverity/Desktop/GRFF_example1.rds")
}



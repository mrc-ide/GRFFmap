
library(tidyverse)
library(rstan)
library(sf)

df_dat <- read.csv("R_ignore/R_scripts/data/Uganda_data.csv") |>
  mutate(t = syear - 2005)

df_dat |>
  ggplot() + theme_bw() +
  geom_point(aes(x = syear, y = x / n))

df_dat |>
  ggplot() + theme_bw() +
  geom_point(aes(x = long, y = lat, col = t), size = 4) +
  scale_color_viridis_c()

df_dat |>
  ggplot() + theme_bw() +
  geom_point(aes(x = long, y = lat, size = (t + 1)^2 / 1, col = x / n)) +
  scale_color_viridis_c()


n_data <- nrow(df_dat)
length_space <- 1
length_time <- 3

# draw features
n_features <- 30
omega_space <- rbind(rnorm(n_features),
                     rnorm(n_features))
omega_time <- rbind(rnorm(n_features))

# make feature maps
X <- df_dat |>
  select(x = long, y = lat) |>
  as.matrix()

feat <- cbind(cos(X %*% omega_space / length_space + df_dat$t %*% omega_time / length_time),
              sin(X %*% omega_space / length_space + df_dat$t %*% omega_time / length_time)) / sqrt(n_features)

# Prepare data list for Stan
data_list <- list(
  N = n_data,
  n_features = n_features,
  feat = feat,
  k = df_dat$x,
  n_trials = df_dat$n,
  psill_shape = 1,
  psill_rate = 1 / 10,
  nugget_shape = 1,
  nugget_rate = 1000,
  mu_mean = 0,
  mu_sd = 10
)

# run stan model
fit <- stan(
  file = "R_ignore/R_scripts/stan/fixedlength_vario.stan",
  data = data_list,
  iter = 1e3,
  chains = 1
)

summary_stats <- summary(fit)$summary
summary_stats[,c("n_eff", "Rhat")]


# Extract samples for beta and std
draws <- extract(fit)
n_draws <- nrow(draws$beta)

hist(draws$psill)
hist(draws$nugget)
hist(draws$mu)

thinning <- round(n_draws / 890)

thin_vec <- seq(1, n_draws, by = thinning)
n_draws_thin <- length(thin_vec)

beta_samples <- draws$beta[thin_vec,]
psill_samples <- draws$psill[thin_vec]
nugget_samples <- draws$nugget[thin_vec]
mu_samples <- draws$mu[thin_vec]




# ------------------------------------------------

# define prediction grid
dx <- 0.1
dy <- 0.1
df_pred <- expand_grid(long = seq(29, 36, dx),
                       lat = seq(-1.5, 4.5, dy))

# read in Uganda shapefile
s <- st_read("/Users/rverity/Dropbox/Bob/Work/Analyses/mapping/GADM/gadm41_UGA_shp/gadm41_UGA_0.shp")

# intersect
coords_sf <- st_as_sf(df_pred, coords = c("long", "lat"), crs = 4326)
df_pred_intersect <- st_join(coords_sf, s, join = st_intersects) |>
  filter(!is.na(COUNTRY))
df_pred <- st_coordinates(df_pred_intersect) |>
  as.data.frame() |>
  rename(long = X, lat = Y)

# make feature maps for prediction
X_pred <- df_pred |>
  select(x = long, y = lat) |>
  as.matrix()

# get fade
t_vec <- seq(0, 16, 0.25)
fade_mat <- get_fade(t_vec, df_dat$t, show_window = 1, fade_window = 1)
#fade_mat <- matrix(1, length(t_vec), n_data)

i <- 60
#for (i in seq_along(t_vec)) {
  message(sprintf("%s of %s", i, length(t_vec)))
  
  # fix prediction at this time
  t_pred <- rep(t_vec[i], nrow(X_pred))
  
  # generate prediction feature maps
  feat_pred <- cbind(cos(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time),
                     sin(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time)) / sqrt(n_features)
  
  #feat_pred <- cbind(cos(X_pred %*% omega_space / length_space),
  #                   sin(X_pred %*% omega_space / length_space)) / sqrt(n_features)
  
  # get predictions over all posterior draws
  post_all <- feat_pred %*% t(beta_samples) |>
    sweep(2, mu_samples, "+") |>
    plogis()
  
  # summarise posterior predictions
  pred_quant <- apply(post_all, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) |>
    t() |>
    as.data.frame()
  
  # get plotting data.frames from infered distribution and truth
  df_pred2 <- df_pred |>
    bind_cols(pred_quant) |>
    mutate(MOE = (`97.5%` - `2.5%`) / 2,
           p_pred = `50%`,
           mask = (MOE > 0.125),
           ex_5 = rowMeans(post_all > 0.05),
           ex_10 = rowMeans(post_all > 0.10))
  
  df_masked <- df_pred2 |>
    filter(mask == TRUE)
  
  plot1 <- df_pred2 |>
    ggplot() + theme_bw() +
    geom_raster(aes(x = long, y = lat, fill = ex_10)) +
    geom_point(aes(x = long, y = lat, fill = x / n, alpha = alpha),
               size = 3, pch = 21, colour = grey(0.5), data = cbind(df_dat, alpha = fade_mat[i,])) +
    #geom_point(aes(x = long, y = lat, fill = x / n), size = 3, pch = 21, colour = grey(0.5), data = cbind(df_dat)) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
    #geom_tile(aes(x = long, y = lat), width = dx, height = dy, fill = grey(0.5), alpha = 0.4, data = df_masked) +
    ggtitle(sprintf("Time = %s", t_vec[i])) +
    coord_fixed(ratio = 1) +
    guides(alpha = "none")
  plot1
  
  # save to file as .png
  ggsave(sprintf("R_ignore/R_scripts/outputs/Uganda/frame%s.png", i),
         plot = plot1, width = 6, height = 5, dpi = 250)
  
#}

# make video
# image_files <- sprintf("R_ignore/R_scripts/outputs/Uganda/frame%s.png", seq_along(t_vec))
# av::av_encode_video(input = image_files,
#                     output = "R_ignore/R_scripts/outputs/Uganda.mp4",
#                     framerate = 20)

library(tidyverse)
library(plotly)
library(rstan)
library(here)

# Define output folder
outFldr <- "R_ignore/R_scripts/outputs/GRFF_africa/"
if (!dir.exists(outFldr)) {
  dir.create(outFldr, recursive = TRUE)
}

sf::sf_use_s2(FALSE)

# Load and process data --------------------------------------------------------
dat <- read.csv("R_ignore/R_scripts/data/validated_and_candidate_get_prevalence.csv")
admin0_africa <- readRDS("R_ignore/R_scripts/data/sf_admin0_africa.rds")  %>%
  sf::st_make_valid() 

union <- admin0_africa %>%
  sf::st_union()

# Filter for individual 675V mutation
dat <- dat %>%
  filter(mutation == "k13:675:V") %>%
  mutate(time = as.Date(sprintf("%s-01-01", year)))

## INTERMEDIATE SOLUTION: aggregate duplicates (or near duplicate) observations
duplicates <- dat %>%
  group_by(latitude, longitude, time) %>%
  group_split()
dat <- dat %>%
  group_by(latitude, longitude, time) %>%
  summarise(
    numerator = sum(numerator),
    denominator = sum(denominator),
    prevalence = numerator / denominator,
    .groups = "drop"
  )

# Transform raw numerator & denominator into z value using empirical logit
dat <- dat %>%
  mutate(z_est = qlogis((numerator + 0.5) / (denominator + 1.0)),
         year = year(time))

# Convert dat into dataframe for model
df_dat <- dat |>
  select(latitude, longitude, time, numerator, denominator, prevalence, z_est) %>%
  mutate(year = year(time)) %>%
  rename(z = z_est,
         p = prevalence,
         N = denominator,
         k = numerator)

# Create grid on which to predict prevalence -----------------------------------
# Create regular grid over Africa
grid <- sf::st_make_grid(union,
                         what = "centers",
                         cellsize = 0.5,   # ~50km at equator
                         square = TRUE)
grid_sf <- sf::st_sf(geometry = grid)

# Filter grid to include only points within Africa
grid <- sf::st_filter(grid_sf, union, .predicate = sf::st_within)

# Extract coordinates (lon, lat)
coords <- sf::st_coordinates(grid) |>
  as_tibble() |>
  rename(lon = X, lat = Y)

# Create space-time grid by crossing with years
df_grid <- expand_grid(coords, year = min(dat$year):2024) |>
  mutate(time = as.Date(sprintf("%s-01-01", year)))

t_vec <- sort(unique(df_grid$time))

# Plot data --------------------------------------------------------------------
# Plot prevalence data
prevalence_admin0_plot <- ggplot() +
  geom_sf(data = admin0_africa, fill = NA, color = "black", linewidth = 0.3) +
  geom_point(data = df_dat %>% arrange(p),
             aes(x = longitude, y = latitude, fill = p),
             shape = 21, size = 1) +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "magma") +
  theme_bw()
prevalence_plot_path <- paste0(outFldr, "prevalence_plot.png")
ggsave(prevalence_plot_path, prevalence_admin0_plot)

# Plot logit transformed data
logit_admin0_plot <- ggplot() +
  geom_sf(data = admin0_africa, fill = NA, color = "black", linewidth = 0.3) +
  geom_point(data = df_dat %>% arrange(p),
             aes(x = longitude, y = latitude, fill = plogis(z)),
             shape = 21, size = 1) +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "magma") +
  theme_bw()
logit_plot_path <- paste0(outFldr, "logit_plot.png")
ggsave(logit_plot_path, logit_admin0_plot)

# Create Spatial and Temporal Feature Maps -------------------------------------
#Parameters
# parameters
length_space <- 2.5
length_time <- 800
sill <- 2
mu <- -2
n_features <- 75

# Draw features
omega_space <- rbind(rnorm(n_features),
                     rnorm(n_features))
omega_time <- rbind(rnorm(n_features))

# Make feature maps
X <- df_dat %>%
  select(longitude, latitude) %>%
  as.matrix()

# Time origin
origin_date <- min(df_dat$time)
t_mat <- as.numeric(difftime(df_dat$time, origin_date, units = "days")) / 365.25
t_mat <- matrix(t_mat, ncol = 1) 

feat <- cbind(cos(X %*% omega_space / length_space + t_mat %*% omega_time / length_time),
              sin(X %*% omega_space / length_space + t_mat %*% omega_time / length_time)) / sqrt(n_features)

# Fit stan model ---------------------------------------------------------------
# prepare data list for Stan
data_list <- list(
  N = dim(df_dat)[1],
  n_features = n_features,
  feat = feat,
  k = df_dat$k,
  n_trials = df_dat$N,
  sill_shape = 1,
  sill_rate = 0.1,
  mu_mean = 0,
  mu_sd = 3
)

# run stan model
fit <- stan(
  file = "R_ignore/R_scripts/stan/fixedlength_v4.stan",
  data = data_list,
  iter = 1e3,
  chains = 1
)

# extract MCMC samples
draws <- extract(fit)
n_draws <- nrow(draws$beta)

beta_samples <- draws$beta
sill_samples <- draws$sill
mu_samples <- draws$mu

hist(sill_samples, breaks = 100)
abline(v = sill, lwd = 2, col = 2)

hist(mu_samples, breaks = 100)
abline(v = mu, lwd = 2, col = 2)

# Obtain predictions --------------------------------------------------------------------

l <- list()
for (i in seq_along(t_vec)) {
  message(sprintf("Predicting for year %s (%s of %s)", year(t_vec[i]), i, length(t_vec)))
  
  # Get subset of df_grid for current time
  df_grid_year <- df_grid %>%
    filter(time == t_vec[i])
  
  # Coordinates
  X_pred <- df_grid_year %>%
    select(lon, lat) %>%
    as.matrix()
  
  # Time feature
  origin_date <- min(dat$time)
  t_pred <- as.numeric(difftime(df_grid_year$time, origin_date, units = "days")) / 365.25
  t_pred <- matrix(t_pred, ncol = 1)  # --> This will now be [79 x 1]
  
  # generate prediction feature maps
  feat_pred <- cbind(cos(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time),
                     sin(X_pred %*% omega_space / length_space + t_pred %*% omega_time / length_time)) / sqrt(n_features)
  
  # get predictions over all posterior draws
  post_all <- feat_pred %*% t(beta_samples) |>
    sweep(MARGIN = 2, STATS = mu_samples, FUN = "+") |>
    plogis()
  
  # summarise posterior predictions
  pred_quant <- apply(post_all, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) |>
    t() |>
    as.data.frame()
  
  # get combined data.frame from inferred distribution and truth
  df_pred <- df_grid_year %>%
    bind_cols(pred_quant) %>%
    mutate(MOE = `97.5%` - `2.5%`,
           p_pred = `50%`,
           mask = (`97.5%` - p_pred) > 0.1 | (p_pred - `2.5%`) > 0.1)
  
  l[[i]] <- df_pred
}
df_combined <- bind_rows(l)

# produce plots
# df_masked <- df_combined |>
#   filter(mask == TRUE) |>
#   mutate(Type = factor("p_pred", levels = c("p_true", "p_pred")))
# dx <- diff(sim_grid$x)[1]
# dy <- diff(sim_grid$y)[1]

ggplot(df_combined) + theme_bw() +
  geom_raster(aes(x = lon, y = lat, fill = p_pred)) +
  geom_point(data = df_dat  %>% arrange(p), aes(x = longitude, y = latitude, fill = p/100),
             shape = 21, color = "grey50", size = 2) +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  ggtitle("Predicted Prevalence by Year for 675V")

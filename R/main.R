
#------------------------------------------------
#' @title Run MCMC for Gaussian Random Feature Field Model
#'
#' @description Converts spatial-temporal data into a random Fourier feature 
#' representation and fits a hierarchical spatial-temporal binomial model 
#' using Stan and MCMC.
#'
#' @param longitude Vector of longitude coordinates (in decimal degrees).
#' @param latitude Vector of latitude coordinates (in decimal degrees).
#' @param survey_day Vector of time points (as Date or numeric).
#' @param n_pos Vector of positive test counts at each site.
#' @param n_samp Vector of total samples at each site.
#' @param length_space Length scale (in meters) for spatial kernel.
#' @param length_time Length scale (in days) for temporal kernel.
#' @param longitude_centre Central longitude for coordinate transformation.
#' @param latitude_centre Central latitude for coordinate transformation.
#' @param time_centre Reference date to center time values.
#' @param mu_mean Prior mean for the global intercept.
#' @param mu_sd Prior standard deviation for the global intercept.
#' @param sigma_shape Shape parameter for prior on GP amplitude.
#' @param sigma_rate Rate parameter for prior on GP amplitude.
#' @param pi_nugget_shape1 Beta prior shape1 for site-level nugget effect.
#' @param pi_nugget_shape2 Beta prior shape2 for site-level nugget effect.
#' @param n_features Number of random Fourier features.
#' @param MCMC_iter Number of MCMC iterations per chain.
#' @param MCMC_chains Number of MCMC chains to run.
#'
#' @return A list containing:
#' \item{data}{Preprocessed data frame with spatial and temporal features.}
#' \item{params}{Model hyperparameters used.}
#' \item{MCMC_fit}{Stan fit object returned by \code{rstan::stan()}.}
#'
#' @export

run_GRFF_mcmc <- function(longitude, latitude, survey_day, n_pos, n_samp,
                          length_space = 75, length_time = 5,
                          longitude_centre = mean(longitude),
                          latitude_centre = mean(latitude),
                          time_centre = min(survey_day),
                          mu_mean = 0, mu_sd = 5,
                          sigma_shape = 1, sigma_rate = 0.1,
                          pi_nugget_shape1 = 1, pi_nugget_shape2 = 9,
                          n_features = 100,
                          MCMC_iter = 1e3,
                          MCMC_chains = 1) {
  
  # Convert longitude and latitude to Cartesian coordinates centered at the specified reference point
  coords_cartesian <- lonlat_to_cartesian(lon = longitude, 
                                          lat = latitude,
                                          lon_centre = longitude_centre, 
                                          lat_centre = latitude_centre)
  
  # Create a dataframe with observed data and transformed coordinates
  df_dat <- data.frame(longitude = longitude,
                       latitude = latitude,
                       survey_day = survey_day,
                       x = coords_cartesian$x,
                       y = coords_cartesian$y,
                       t = as.numeric(survey_day - time_centre),  # days since time_centre
                       n_pos = n_pos,
                       n_samp = n_samp)
  
  # Draw random Fourier frequencies for space and time
  omega_space <- rbind(rnorm(n_features),
                       rnorm(n_features))  # 2 x n_features for 2D spatial features
  omega_time <- rbind(rnorm(n_features))    # 1 x n_features for time
  
  # Prepare input matrices for the feature map
  X_space <- cbind(df_dat$x, df_dat$y)
  X_time <- df_dat$t
  
  # Construct random Fourier feature representation (cosine and sine basis)
  feat <- cbind(cos(X_space %*% omega_space / length_space + X_time %*% omega_time / length_time),
                sin(X_space %*% omega_space / length_space + X_time %*% omega_time / length_time)) / sqrt(n_features)
  
  # Assemble data list to pass to Stan model
  data_list <- list(
    n_sites = nrow(df_dat),
    n_features = n_features,
    feat = feat,
    n_pos = df_dat$n_pos,
    n_samp = df_dat$n_samp,
    mu_mean = mu_mean,
    mu_sd = mu_sd,
    sigma_shape = sigma_shape,
    sigma_rate = sigma_rate,
    pi_nugget_shape1 = pi_nugget_shape1,
    pi_nugget_shape2 = pi_nugget_shape2
  )
  
  # Compile and run the Stan model using the preprocessed data
  fit <- stan(
    file = system.file("stan", "RFF_v1.stan", package = "GRFFmap"),
    data = data_list,
    iter = MCMC_iter,
    chains = MCMC_chains
  )
  
  # Package model hyperparameters into a list for return
  params <- list(mu_mean = mu_mean,
                 mu_sd = mu_sd,
                 sigma_shape = sigma_shape,
                 sigma_rate = sigma_rate,
                 pi_nugget_shape1 = pi_nugget_shape1,
                 pi_nugget_shape2 = pi_nugget_shape2,
                 n_features = n_features,
                 MCMC_iter = MCMC_iter,
                 MCMC_chains = MCMC_chains)
  
  # Return the full model output, parameters, and processed data
  ret <- list(data = df_dat,
              params = params,
              MCMC_fit = fit)
  return(ret)
}

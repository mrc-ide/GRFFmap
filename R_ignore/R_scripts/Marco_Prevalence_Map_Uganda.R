
library(mapview)
library(INLA)
library(tidyverse)
library(cccd)

df <- read_csv("R_ignore/R_scripts/data/Uganda_data_cleaned.csv")

length(which(is.na(df$prev) == FALSE)) # total number of observed prevalence points in UGA from 2010-2021

# time variable
df$time = rep(1:n_distinct(df$syear),
              n_distinct(df$sitesID))

# Interactive map
sites = unique(df[,c('sitesID','long', 'lat')])
mapview(sites, xcol = "long", ycol = "lat", crs = 4326, grid = FALSE)

# mapping
ggplot(df) + theme_bw() +
  geom_point(aes(x = long, y = lat, color = prev, size = n)) +
  coord_fixed(ratio = 1) +
  scale_color_gradient(low = "blue", high = "orange")

# 
coords = unique(df[c("long", "lat")]) # coordinates for each site
boundary = inla.nonconvex.hull(as.matrix(coords[,1:2]))
mesh = inla.mesh.2d(boundary = boundary, max.edge = c(0.8, 1.3), cutoff = 0.1)
ggplot() +
  gg(mesh) +
  geom_point(data = coords, aes(long, lat))

# 
# approach 1: straigtforward
# 
# create SPDE without specifying detailed priors
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
spde$n.spde

# create an index set for the SPDE model to help INLA interpret the vertices of the mesh correctly.
n_years = length(unique(df$time))
n_spatial = mesh$n
s_index = inla.spde.make.index(
  name = "spatial.field",
  n.spde = n_spatial,
  n.group = n_years)

# create a projection matrix to project the GRF from the observations to the vertices of our mesh
# for this we need the coordinates of our data in matrix form
coordinates.allyears = df[c("long", "lat")] %>% as.matrix()
A_est = inla.spde.make.A(mesh = mesh,
                         loc = coordinates.allyears,
                         group = df$time,
                         n.group = n_years)


dim(A_est) # 13419 82845
length(s_index$spatial.field) # 82845

# stack for estimation
stack_est <- inla.stack(
  tag = "est",
  data = list(Npos = df$x,
              Ntrials = df$n),
  A = list(1, A_est),
  effects = list(data.frame(b0 = rep(1, nrow(df))), s = s_index)
)

# import shapefile of Uganda
AFR <- st_read("/Users/rverity/Dropbox/Bob/Work/Analyses/mapping/GADM/gadm41_UGA_shp/gadm41_UGA_0.shp")
#AFR <- st_read("~/Downloads/result/8e76854e-6afa-4b3e-9efa-42594d22176c/afr_g2014_2013_0.shp")
uganda_shp <- AFR
#uganda_shp <- AFR %>% filter(ISO3 == "UGA")
#uganda_shp <- st_transform(uganda_shp, crs = "+proj=longlat +datum=WGS84")

# Create a grid of points over the bounding box of Uganda
uganda_bbox <- st_bbox(uganda_shp)
long_range <- seq(uganda_bbox["xmin"], uganda_bbox["xmax"], by = 0.05) # 0.15 degrees = 18km
lat_range <- seq(uganda_bbox["ymin"], uganda_bbox["ymax"], by = 0.05)
prediction_times <- unique(df$time)

# Create a grid of points
grid_points <- expand.grid(long = long_range, lat = lat_range, time = prediction_times)

# Convert grid points to an sf object
grid_sf <- st_as_sf(grid_points, coords = c("long", "lat"), crs = "+proj=longlat +datum=WGS84")

# Clip the grid to the boundary of Uganda
grid_sf <- st_intersection(grid_sf, uganda_shp)

# convert grid_sf back to long lat
prediction_df <- as.data.frame(st_coordinates(grid_sf))
prediction_df$time <- grid_sf$time

# visualise the prediction grid
# prediction_df %>%
#   distinct(X, Y) %>%
#   ggplot()+
#   geom_point(aes(X, Y)) +
#   gg(mesh)

# Create the projection matrix for prediction
coordinates_pred <- prediction_df[c("X", "Y")] %>% as.matrix()
A_pred <- inla.spde.make.A(
  mesh = mesh,
  loc = coordinates_pred,
  group = prediction_df$time,
  n.group = n_years
)

# Create the stack for prediction
stack_pred <- inla.stack(
  tag = "pred",
  data = list(Npos = NA, Ntrials = NA),
  A = list(1, A_pred),
  effects = list(data.frame(b0 = rep(1, nrow(prediction_df))), s = s_index)
)

# Combine the estimation and prediction stacks
stack_full <- inla.stack(stack_est, stack_pred)

# add time component to the model using ar1
# AR1 prior
rho_hyper = list(theta = list(prior = "pccor1", param = c(0.9, 0.4)))

# INLA formula
f = Npos ~ 0 + b0 + 
  f(spatial.field, 
    model = spde, 
    group = spatial.field.group, 
    control.group = list(model = "ar1", hyper = rho_hyper))

# Run INLA with the updated formula
t0 <- Sys.time()
results = inla(f, family = "binomial", Ntrials = Ntrials,
               control.family = list(link = "logit"),
               data = inla.stack.data(stack_full),
               control.predictor = list(
                 compute = TRUE, link = 1,
                 A = inla.stack.A(stack_full)),
               control.compute = list(
                 config = TRUE,
                 return.marginals.predictor = TRUE)
)
Sys.time() - t0

summary(results)

# Extract and visualize the predictions
index_pred <- inla.stack.index(stack_full, "pred")$data
predictions <- data.frame(
  long = prediction_df$X,
  lat = prediction_df$Y,
  time = prediction_df$time,
  mean = results$summary.fitted.values[index_pred, "mean"],
  sd = results$summary.fitted.values[index_pred, "sd"]
)

# Visualize the predictions
ggplot(predictions, aes(x = long, y = lat, fill = mean)) +
  geom_tile() +
  scale_fill_viridis() +
  ggtitle("Predicted Prevalence") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  facet_wrap(~ time, ncol = 3) +
  geom_point(aes(x = long, y = lat, fill = x / n), pch = 21, data = subset(df, !is.na(df$n)))




# extracting the posterior samples
names(results$marginals.fixed)
names(results$marginals.random)
names(results$marginals.hyperpar)
round(results$summary.fixed[,c("mean","0.025quant","0.975quant")],3)
modfix = results$summary.fixed
modfix

# density plot for model intercept
plot(results$marginals.fix[[1]],type ='l',xlab=expression(beta[0]),ylab="density")
abline(v = modfix[1, c(3, 5)], lty=2)
  # The model shows a low estimated log-odds of the baseline ART-R prevalence across all time points = -3.715 (95% CrI: -4.752 to -2.794)
  # odds = exp(-3.715)
  # prob= odds / (odds+1) = 0.024 = 2.4% (prob of ART-R)

# posterior estimates of hyperparam
modhy = results$summary.hyperpar
modhy
output.field = inla.spde2.result(inla = results,
                                 name = "spatial.field",
                                 spde = spde,
                                 do.transf = TRUE)
out.range = exp(output.field$summary.log.range.nominal)
out.var = exp(output.field$summary.log.variance.nominal)

#par(mfrow=c(2,2))

# AR1 parameter (magnitude of temporal correlation)
plot(results$marginals.hyperpar$`GroupRho for spatial.field`,type = 'l',xlab = expression(rho),ylab = "density", xlim=c(0.6,1))

# variance (how much spatial variation) very uncertain in the spatial field
plot(output.field$marginals.variance.nominal[[1]],type = 'l',xlab = expression(sigma^2),ylab = "density")

# range (how far-reaching the spatial dependencies are) very strong up to a distance of 2-3 degrees of long/lat
plot(output.field$marginals.range.nominal[[1]],type = 'l',xlab = "spatial range",ylab = "density")

# compare predicted vs observed
# comparing predicted with observed
hist(predictions$mean, breaks = 20)
hist(df$prev, breaks = 20)

idat = inla.stack.index(stack_est, 'est')$data
cor(df$prev, exp(results$summary.linear.predictor$mean[idat]), use="complete.obs")
plot(df$prev, exp(results$summary.linear.predictor$mean[idat]), xlab="predicted", ylab="observed")












# exceedance probability
EP_thresh <- 0.1
EP_prob_thresh <- 0.8

t0 <- Sys.time()
exceedance_prob <- mapply(function(x) {
  1 - inla.pmarginal(q = EP_thresh, marginal = x)
}, results$marginals.fitted.values[index_pred])
Sys.time() - t0

predictions$EP <- exceedance_prob


ggplot(predictions, aes(x = long, y = lat, fill = EP)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", name = "Probability\nprevalence > 10%") +
  coord_fixed(ratio = 1) +
  facet_wrap(~ time, ncol = 3) +
  theme_minimal() +
  ggtitle("Exceedance Probability of Prevalence > 10%")

high_prob_areas <- exceedance_prob >= EP_prob_thresh
predictions$high_prob_80 <- high_prob_areas

magma_colors <- viridis::viridis(2, option = "viridis")

ggplot(predictions, aes(x = long, y = lat, fill = high_prob_80)) +
  geom_tile() +
  scale_fill_manual(values = magma_colors, name = "High Probability\n(> 80%)") +
  coord_fixed(ratio = 1) +
  facet_wrap(~ time, ncol = 3) +
  theme_minimal() +
  ggtitle("High Exceedance Probability of Prevalence > 10% (Threshold: 80%)")

# areas of interest
interest <- predictions[predictions$time == 8 & predictions$high_prob_80 == TRUE,]
# identify districts of areas of interest
data(world.cities)
  # Convert the points of interest to an sf object
interest_sf <- st_as_sf(interest, coords = c("long", "lat"), crs = 4326)
  # Convert the world.cities dataset to an sf object
cities_sf <- st_as_sf(world.cities, coords = c("long", "lat"), crs = 4326)
  # find the nearest city for each point
nearest_city_indices <- st_nearest_feature(interest_sf, cities_sf)
nearest_cities <- cities_sf[nearest_city_indices, ]
interest$city <- nearest_cities$name
  # Plot
plot(st_geometry(uganda_shp), col = "lightgrey") + plot(st_geometry(interest_sf), add = TRUE, col = "red", pch = 16)
plot_a <- ggplot() +
  geom_sf(data = uganda_shp, fill = "lightgrey", color = NA) +
  geom_sf(data = interest_sf, color = "red", size = 3) +
  theme_minimal() +
  labs(title = "Kitgum, Uganda")
  # find nearest point to Kitgum
kitgum_coords <- world.cities[world.cities$name == "Kitgum", c("long", "lat")]
distances <- sqrt((interest$long - kitgum_coords$long)^2 + (interest$lat - kitgum_coords$lat)^2)
kitgum_index <- which.min(distances)
rownames(interest)[kitgum_index]
  # identify row in predictions dataframe closest to Kitgum
predictions[as.numeric(rownames(interest)[kitgum_index]),]
  # test
ye <- predictions[predictions$long == predictions[as.numeric(rownames(interest)[kitgum_index]),]$long & 
                    predictions$lat == predictions[as.numeric(rownames(interest)[kitgum_index]),]$lat,]
testing_1 <- results$marginals.fitted.values[[inla.stack.index(stack_full, "pred")$data[as.numeric(rownames(ye))[1]]]]
testing_6 <- results$marginals.fitted.values[[inla.stack.index(stack_full, "pred")$data[as.numeric(rownames(ye))[6]]]]
testing_8 <- results$marginals.fitted.values[[inla.stack.index(stack_full, "pred")$data[as.numeric(rownames(interest)[kitgum_index])]]]
testing_12 <- results$marginals.fitted.values[[inla.stack.index(stack_full, "pred")$data[as.numeric(rownames(ye))[12]]]]
plot_1 <- ggplot(testing_1) + geom_density(aes(x = x))+ xlim(0, 1)
plot_6 <- ggplot(testing_6) + geom_density(aes(x = x)) + xlim(0, 1)
plot_8 <- ggplot(testing_8) + geom_density(aes(x = x)) + xlim(0, 1)
plot_12 <- ggplot(testing_12) + geom_density(aes(x = x)) + xlim(0, 1)
gridExtra::grid.arrange(plot_a,plot_1, plot_6, plot_8, plot_12, ncol = 1)
# ribbon plot?
# sampling sites plot


# Cross validation - testing

# training vs validation sets
set.seed(23876)
validation_sites <- sample(unique(df$sitesID), size = 5)  # 5 sites for validation
training_data <- df[!df$sitesID %in% validation_sites, ]
validation_data <- df[df$sitesID %in% validation_sites, ]

# mesh (training)
coords_train <- unique(training_data[c("long", "lat")])
boundary_train <- inla.nonconvex.hull(as.matrix(coords_train))
mesh_train <- inla.mesh.2d(boundary = boundary_train, max.edge = c(0.8, 1.3), cutoff = 0.1)

# SPDE model (training)
spde_train <- inla.spde2.matern(mesh = mesh_train, alpha = 2, constr = TRUE)

# index set (training)
n_years <- length(unique(training_data$time))
s_index_train <- inla.spde.make.index(name = "spatial.field", n.spde = mesh_train$n, n.group = n_years)

# projection matrix (training)
A_est_train <- inla.spde.make.A(mesh = mesh_train, loc = as.matrix(training_data[c("long", "lat")]), group = training_data$time, n.group = n_years)

# Estimation stack (training)
stack_est_train <- inla.stack(
  tag = "est",
  data = list(Npos = training_data$x, Ntrials = training_data$n),
  A = list(1, A_est_train),
  effects = list(data.frame(b0 = rep(1, nrow(training_data))), s = s_index_train)
)

# prediction stack (validation)
A_pred_val <- inla.spde.make.A(mesh = mesh_train, 
                               loc = as.matrix(validation_data[c("long", "lat")]), 
                               group = validation_data$time, 
                               n.group = n_years)

stack_pred_val <- inla.stack(
  tag = "pred",
  data = list(Npos = NA, Ntrials = NA),
  A = list(1, A_pred_val),
  effects = list(data.frame(b0 = rep(1, nrow(validation_data))), s = s_index_train)
)

# combine the stacks
stack_full <- inla.stack(stack_est_train, stack_pred_val)

# define AR1 prior
rho_hyper = list(theta = list(prior = "pccor1", param = c(0.9, 0.4)))

# formula
f <- Npos ~ 0 + b0 + 
  f(spatial.field, model = spde_train, group = spatial.field.group, control.group = list(model = "ar1", hyper = rho_hyper))

# run INLA
t0 <- Sys.time()
results_cv <- inla(f, family = "binomial", Ntrials = Ntrials,
                control.family = list(link = "logit"),
                data = inla.stack.data(stack_full),
                control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stack_full))
                )
Sys.time() - t0

# Compare predicted vs observed
index_pred <- inla.stack.index(stack_full, "pred")$data
predictions <- data.frame(
  long = validation_data$long,
  lat = validation_data$lat,
  time = validation_data$time,
  observed = validation_data$prev,
  predicted = results$summary.fitted.values[index_pred, "mean"],
  sd = results$summary.fitted.values[index_pred, "sd"]
)
predictions <- predictions[!is.na(predictions$observed),]

# Plot observed vs predicted values
ggplot(predictions) +
  geom_point(aes(x = observed, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  ggtitle("Observed vs Predicted Prevalence")

# Calculate correlation
cor(predictions$observed, predictions$predicted, use = "complete.obs")













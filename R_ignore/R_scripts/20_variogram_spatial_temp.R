# Packages
library(sf)
library(sp)
library(spacetime)
library(gstat)
library(dplyr)
library(lubridate)

# read in prevalence data
dat <- read.csv("R_ignore/R_scripts/data/all_who_get_prevalence_africa.csv") |>
  mutate(collection_day = as.Date(collection_day),
         collection_year = year(collection_day)) |>
  select(study_id, study_name, site_name, latitude, longitude, collection_year, prevalence)


# Clamp prevalence strictly inside (0,1) to avoid Inf on logit
eps <- 1e-6
dat <- dat %>%
  mutate(
    prevalence = pmin(pmax(prevalence, eps), 1 - eps),
    date       = as.Date(paste0(collection_year, "-07-01")),
    prev_logit = log(prevalence / (1 - prevalence))
  )

# ---- Build a spacetime object (STFDF) ----
dat2 <- dat %>%
  rename(long = longitude, lat = latitude) %>%
  mutate(time = as.POSIXct(date, tz = "UTC"))

coords  <- dat2 %>% select(long, lat)
sp      <- SpatialPointsDataFrame(coords,
                                  data = data.frame(region = dat2$site_name),
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
time    <- dat2$time
values  <- data.frame(prevalence = dat2$prev_logit)

# Sanity checks (must all be equal)
stopifnot(nrow(sp) == nrow(values), length(time) == nrow(values))

# Create STIDF
stidf <- STIDF(sp = sp, time = time, data = values)

df <- as(stidf, "data.frame")  # columns: time, x, y, plus attributes
df <- df[is.finite(df$prevalence) & !is.na(df$time), c("time","x","y","prevalence")]

sp_clean <- SpatialPointsDataFrame(df[, c("x","y")],
                                   data = data.frame(id = seq_len(nrow(df))),
                                   proj4string = CRS(proj4string(stidf@sp)))

stidf_ok <- STIDF(sp = sp_clean,
                  time = as.POSIXct(df$time, origin="1970-01-01", tz="UTC"),
                  data = data.frame(prevalence = df$prevalence))

# Choose spatial binning (example: up to 300 km, in 25 km bins)
cutoff_m <- 300e3
width_m  <- 25e3

# 0 to 5-year lags, expressed in days
yr  <- 365.25
tl  <- 0:5 * yr

v_emp <- variogramST(
  prevalence ~ 1, data = stidf,
  tlags = tl,
  tunit = "days",          # << key: interpret tlags as days
  cutoff = 300e3, width = 25e3,
  na.omit = TRUE
)

plot(vgm_st, map = FALSE)

# ---- 3) Sample spatio-temporal variogram ----
# Choose binning for space (meters) and time (days)
vgm_st <- variogramST(
  prev_logit ~ 1, data = st_obj,
  tlags = tlags_days,
  tunit = "days",             # ensures tlags are read in days
  assumeRegular = FALSE,
  cutoff = 500000,            # spatial cutoff in meters
  width  = 50000              # spatial bin width in meters
)

# Quick look
plot(vgm_st, map = FALSE)

# ---- 4) Fit a simple spatio-temporal variogram model ----
# Start from marginal (pure) spatial + temporal models
vgm_S <- vgm(psill = 0.8, model = "Exp", range = 200000, nugget = 0.2)  # spatial
vgm_T <- vgm(psill = 0.8, model = "Exp", range = 365*1.0, nugget = 0.0) # temporal (~1 year)

# Separable model (common, interpretable)
sep_model <- vgmST("separable",
                   space = vgm_S,
                   time  = vgm_T,
                   sill  = 1.0)  # overall sill scaling

fit_sep <- fit.StVariogram(vgm_st, sep_model, method = "L-BFGS-B", lower = 1e-6)
fit_sep

# Plot fitted vs empirical
plot(vgm_st, fit_sep, wireframe = TRUE)

# ---- 5) (Optional) Use in spatio-temporal kriging ----
# newdata_st: STF for prediction grid x times; then:
# kr <- krigeST(prev_logit ~ 1, data = st_obj, newdata = newdata_st, modelList = fit_sep)
# Back-transform predictions if you used logit: inv_logit <- function(x) 1/(1+exp(-x))

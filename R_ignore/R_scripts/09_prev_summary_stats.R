#---------------Summary Stats for Modeling Output-------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(here)
  library(sf)
  library(devtools)
  library(dplyr)
  library(fs)
  library(here)
  library(lwgeom)
})

load_all()

set.seed(1)
sf_use_s2(FALSE)

# --------------------------- Settings --------------------------------------
# prediction parameters
nx <- 200
ny <- 200
t_vec <- 1995:2025
t_num <- length(t_vec)

# plotting parameters
plot_times <- seq(2010, 2024, by = 1) # should be within t_vec
t_win <- 2           # window around each facet time (years). Points within this window are displayed


# --------------------------- Load & filter data ----------------------------

CACHE_DIR <- "R_ignore/R_scripts/outputs/model_outputs/supplemental/GRFF_model_output_all_WHO_mutations/"
dir_create(CACHE_DIR)

# read in prevalence data
dat <- read.csv("R_ignore/R_scripts/data/all_mutations_get_prevalence.csv") |>
  mutate(collection_day = as.Date(collection_day)) |>
  select(study_id, survey_id, longitude, latitude, year, country_name, numerator, denominator, prevalence, mutation)

# read Africa shape files
africa_shp_admin1 <- readRDS(file = "R_ignore/R_scripts/data/sf_admin1_africa.rds")


# --------------------------- Fliter admin regions ----------------------------
target_crs <- 4326
africa_shp_admin1 <- st_transform(africa_shp_admin1, crs = target_crs)

# Define East Africa box and crop Admin1 shapefile
bbox_east_africa <- get_east_africa_bbox(target_crs)

# Using the cropped Admin 1 data for East Africa
admin1_regions <- st_make_valid(africa_shp_admin1) |> st_crop(bbox_east_africa)

# --- Add K13 overall data -----------------------------------------------------
dat_with_k13 <- add_combined_k13(dat)

# --------------------------- Define K13 mutants ----------------------------
# which mutation we are focusing on
all_who_mutations <- c("k13:comb", "k13:675:V", "k13:622:I", "k13:469:Y", "k13:446:I", "k13:458:Y", "k13:476:I",   "k13:493:H",   "k13:539:T",
                       "k13:543:T",  "k13:553:L",   "k13:561:H",   "k13:574:L",  "k13:580:Y",
                       "k13:441:L", "k13:449:A",   "k13:469:F",   "k13:481:V",
                       "k13:515:K", "k13:527:H",  "k13:537:I", "k13:537:D", "k13:538:V",  "k13:568:G")
#all_who_mutations = c("k13:675:V", "k13:622:I", "k13:561:H")

all_sum_stats <- list()
for (mut in all_who_mutations){
  cached <- make_or_load(mut)
  p_post_median <- cached$p_post_median
  
  xs <- cached$xs
  ys <- cached$ys
  
  plot_idx <- match(plot_times, t_vec)  # indices in 1:t_num
  p_long_median <- make_raster_long(p_post_median[ , , plot_idx, drop = FALSE], xs, ys, plot_times)
  
  p_median_sf <- st_as_sf(p_long_median, coords = c("x", "y"), crs = 4326, remove = FALSE)
  
  p_median_admin1 <- st_join(p_median_sf, africa_shp_admin1) %>% 
    select(!c(type_2, type_3, id_2, id_3, name_2, name_3))
  
  avg_admin1_year <- p_median_admin1 %>%
    st_drop_geometry() %>%                 # drop point geometry for summarise
    filter(!is.na(id_1)) %>%               # keep only points that matched an admin1
    group_by(id_0, id_1, name_1, name_0, t) %>%     # t is your year/time label from plot_times
    summarise(
      mean_p = mean(p*100, na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(mut = mut)
  
  overall_prev_per_year <- avg_admin1_year %>%
    filter(t %in% c(2012, 2024)) %>%
    pivot_wider(names_from = t, values_from = mean_p, names_prefix = "y_") %>%
    mutate(prev_per_year = (y_2024 - y_2012) / (2024 - 2012) *100) %>%
    select(id_0, id_1, name_1, name_0, prev_per_year)
  
  #all_sum_stats <- all_sum_stats %>% bind_rows(results)
  
  write.csv(avg_admin1_year, paste0("R_ignore/R_scripts/outputs/summary_stats/",mut,"_avg_admin1_year.csv"),row.names = FALSE)
  write.csv(overall_prev_per_year, paste0("R_ignore/R_scripts/outputs/summary_stats/",mut,"_overall_prev_per_year.csv"),row.names = FALSE)
  
}


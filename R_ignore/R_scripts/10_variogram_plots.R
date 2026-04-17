# ============================================================
# Plotting variogram 
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

load_all()

# --------------------------- Deinfe file paths ----------------------------
VARIOGRAM_DIR <- "R_ignore/R_scripts/outputs/model_outputs/variogram_distances"
CACHE_DIR <- "supplemental/variogram"

mutations <- c("k13:comb", mut = "crt76T_mdr86Y")

for (mut in mutations){
  # Spatial variogram
  sp_dat <- readRDS(file.path(VARIOGRAM_DIR, paste0(mut, "_spatial_vg_plot_data.rds")))
  ell_val <- format(round(sp_dat$ell_km, 1), nsmall = 1)
  spatial_variogram <- ggplot() +
    geom_point(data = sp_dat$obs_df, aes(x = dist, y = gamma),
               shape = 1, colour = "black") +
    geom_line(data = sp_dat$vg_df, aes(x = dist, y = fitted), linewidth = 0.8) +
    labs(
      x = "Distance (km)",
      y = "Semivariance",
      title = bquote("Spatial variogram, " ~ ell[km] == .(ell_val))
    ) +
    theme_bw()
  
  # Save figures
  save_figs(file.path(CACHE_DIR, paste0(mut, "_variogram_spatial")), spatial_variogram)
}


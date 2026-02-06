#-----------------------------------------------------------------------------
# Code 20: Trait Divergence Mapping on Environmental Space (Fig 4b)
#
# Description:
#   Visualizes the Standardized Effect Size (SES) of Taxonomic Mean Pairwise 
#   Distance (TaxMPD) within the Environmental Space.
#   Uses GAM (Generalized Additive Models) to smooth surface contours.
#
# Inputs:
#   1. Global_Ecoregion_Environment_EnvPCs.csv
#   2. Diversity_Taxonomic_Traits_SES_Summary.csv (SES Results from Code 19)
#   3. PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv (Reference ID list)
#
# Outputs:
#   - Fig4_b.pdf
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

# [Visualization Settings]
GAM_K       <- 15           # Smoothing knots (Surface complexity)
GRID_RES    <- 300          # Resolution of the raster grid
CUSTOM_COLS <- c("#427638", "#9DC86E", "#FDFBE9", "#D27EAC", "#8C1C54")

# Load Packages
pkgs <- c("data.table", "ggplot2", "mgcv", "sp", "scales")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(ggplot2)
library(mgcv)       
library(sp)         
library(scales)     

cat(">>> 1. Loading and Preparing Data...\n")

# 2. Data Loading
#-----------------------------------------------------------------------------

# A. Load Pre-calculated Environmental PCs
#    (Contains ECO_ID, EnvPC1, EnvPC2...)
env_pc_path <- "Global_Ecoregion_Environment_EnvPCs.csv"
if (!file.exists(env_pc_path)) stop("Error: 'Global_Ecoregion_Environment_EnvPCs.csv' not found.")
df_env_pcs <- fread(env_pc_path)

# B. Load Trait Diversity SES Results
trait_path <- "Diversity_Taxonomic_Traits_SES_Summary.csv"
if (!file.exists(trait_path)) stop("Error: 'Diversity_Taxonomic_Traits_SES_Summary.csv' not found.")
df_trait <- fread(trait_path)

# C. Load Reference File (for Filtering Ecoregions)
#    Updated filename as requested
ref_path <- "PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv"
if (!file.exists(ref_path)) {
  stop(paste("Error: Reference file", ref_path, "not found!"))
}
df_ref <- fread(ref_path)
valid_ids_ref <- unique(df_ref$ECO_ID)

cat(sprintf("   -> Reference IDs loaded: %d\n", length(valid_ids_ref)))

# 3. Merging and Filtering
#-----------------------------------------------------------------------------

# Select relevant columns from SES results (Only TaxMPD is needed for Fig 4b)
trait_subset <- df_trait[, .(ECO_ID, TaxMPD_SES)]

# Merge EnvPCs with SES Data
merged_data <- merge(df_env_pcs, trait_subset, by = "ECO_ID")

# Remove NAs
merged_data <- na.omit(merged_data)
count_before <- nrow(merged_data)

# [CRITICAL] Filter to match the reference Ecoregions
merged_data <- merged_data[ECO_ID %in% valid_ids_ref]
count_after <- nrow(merged_data)

cat(sprintf("   -> Data Filtering:\n"))
cat(sprintf("      - Before filtering: %d\n", count_before))
cat(sprintf("      - After matching reference: %d ( Analysis Set)\n", count_after))

if (count_after == 0) stop("Error: No common ECO_IDs found!")

# Ensure EnvPC1 and EnvPC2 are present
if (!all(c("EnvPC1", "EnvPC2") %in% names(merged_data))) {
  stop("Error: EnvPC1 or EnvPC2 missing from environmental file.")
}

# 4. Generate GAM Surface Plot for TaxMPD SES
#-----------------------------------------------------------------------------
cat(">>> 2. Generating TaxMPD SES Plot...\n")

# A. Prepare Fitting Data
df_fit <- merged_data[, .(EnvPC1, EnvPC2, Trait_Value = TaxMPD_SES)]

# B. GAM Fitting (Value ~ s(PC1, PC2))
#    Using a tensor product smooth or thin plate regression spline
gam_model <- gam(Trait_Value ~ s(EnvPC1, EnvPC2, k = GAM_K), data = df_fit)

# Extract R-squared
r2_val <- summary(gam_model)$r.sq
title_str <- sprintf("TaxMPD (SES) (Adj. R² = %.3f)", r2_val)
cat(sprintf("   -> Model Fitted. R² = %.3f\n", r2_val))

# C. Create Prediction Grid
x_rng <- range(df_fit$EnvPC1)
y_rng <- range(df_fit$EnvPC2)
grid_data <- expand.grid(
  EnvPC1 = seq(x_rng[1], x_rng[2], length.out = GRID_RES),
  EnvPC2 = seq(y_rng[1], y_rng[2], length.out = GRID_RES)
)

# Predict values
grid_data$Pred_Value <- predict(gam_model, newdata = grid_data)

# D. Masking (Convex Hull)
#    Prevent plotting outside the data range (interpolation only, no extrapolation)
hull_idx <- chull(df_fit$EnvPC1, df_fit$EnvPC2)
hull_pts <- df_fit[hull_idx, c("EnvPC1", "EnvPC2")]

# Check if grid points are inside the hull
in_poly <- point.in.polygon(grid_data$EnvPC1, grid_data$EnvPC2, 
                            hull_pts$EnvPC1, hull_pts$EnvPC2)
grid_data$Pred_Value[in_poly == 0] <- NA

# E. Color Scaling
val_min <- min(grid_data$Pred_Value, na.rm=TRUE)
val_max <- max(grid_data$Pred_Value, na.rm=TRUE)

# Rescale for gradient mapping (Min -> Mid-Low -> 0 -> Mid-High -> Max)
rescale_vals <- scales::rescale(c(val_min, val_min/2, 0, val_max/2, val_max))

# F. Plotting
p_taxmpd <- ggplot() +
  # 1. Raster Layer (Smooth Surface)
  geom_raster(data = grid_data, aes(x = EnvPC1, y = EnvPC2, fill = Pred_Value), interpolate = TRUE) +
  
  # 2. Contour Lines
  geom_contour(data = grid_data, aes(x = EnvPC1, y = EnvPC2, z = Pred_Value), 
               colour = "grey50", linewidth = 0.3, bins = 6, na.rm = TRUE) +
  
  # 3. Color Scale
  scale_fill_gradientn(
    colors = CUSTOM_COLS,
    values = rescale_vals,
    name = "SES Value",
    na.value = "transparent"
  ) +
  
  # 4. Hull Boundary (White Border)
  geom_polygon(data = hull_pts, aes(x = EnvPC1, y = EnvPC2), 
               fill = NA, colour = "white", linewidth = 1.2) +
  
  # 5. Theme and Layout
  labs(title = title_str, x = "EnvPC1", y = "EnvPC2") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    aspect.ratio = 1, # Force Square Plot
    
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title = element_text(color = "black", size = 12),
    axis.text  = element_text(color = "black", size = 10),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    
    legend.position = "right",
    legend.key.height = unit(1.0, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.title = element_text(size = 10)
  )

# 5. Save Output
#-----------------------------------------------------------------------------
out_file <- "Fig4_b.pdf"

# Save as PDF (Square dimensions for single plot)
ggsave(out_file, p_taxmpd, width = 7, height = 7, dpi = 300)

cat(sprintf("\n>>> Success! Fig4_b saved to: %s\n", file.path(getwd(), out_file)))

#-----------------------------------------------------------------------------
# Code 21: Trait Divergence Mapping on Trait Space (Fig 4d)
#
# Description:
#   Visualizes the Standardized Effect Size (SES) of Taxonomic Mean Pairwise 
#   Distance (TaxMPD) within the Global Trait Space (PC1 vs PC2).
#   Uses GAM (Generalized Additive Models) to smooth surface contours.
#
# Inputs:
#   1. PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv (Trait Space Coords)
#   2. Diversity_Taxonomic_Traits_SES_Summary.csv (SES Results)
#
# Outputs:
#   - Fig4_d.pdf
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

# [Visualization Settings]
GAM_K       <- 15           # Smoothing knots
GRID_RES    <- 300          # Grid resolution (300x300)

# [Color Palette]
# Blue (Low) -> White (Mid) -> Pink (High) - Consistent with Trait Space visuals
CUSTOM_COLS <- c("#3781B0", "white", "#D27EAC")

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

# A. Load Trait Space Coordinates (PC1, PC2)
coord_path <- "PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv"
if (!file.exists(coord_path)) stop(paste("Error:", coord_path, "not found."))
df_coords <- fread(coord_path)

# B. Load SES Metrics
ses_path <- "Diversity_Taxonomic_Traits_SES_Summary.csv"
if (!file.exists(ses_path)) stop(paste("Error:", ses_path, "not found."))
df_ses <- fread(ses_path)

# 3. Merging and Processing
#-----------------------------------------------------------------------------

# A. Extract Coordinates
coord_subset <- df_coords[, .(ECO_ID, PC1, PC2)]

# B. Extract SES Metrics (Only TaxMPD is required for Fig 4d)
ses_subset <- df_ses[, .(ECO_ID, TaxMPD_SES)]

# C. Merge
df_analysis <- merge(coord_subset, ses_subset, by = "ECO_ID")
df_analysis <- na.omit(df_analysis)

cat(sprintf("   -> Analysis based on %d ecoregions.\n", nrow(df_analysis)))

# 4. Generate GAM Surface Plot for TaxMPD SES
#-----------------------------------------------------------------------------
cat(">>> 2. Generating TaxMPD SES Plot (Trait Space)...\n")

# A. Prepare Data
# Using PC1 and PC2 (Trait Space)
df_fit <- df_analysis[, .(PC1, PC2, Trait_Value = TaxMPD_SES)]

# B. GAM Fitting
#    Model: Value ~ Smooth(PC1, PC2)
gam_model <- gam(Trait_Value ~ s(PC1, PC2, k = GAM_K), data = df_fit)

# Extract R-squared
r2_val <- summary(gam_model)$r.sq
title_str <- sprintf("TaxMPD (SES) (Adj. R² = %.3f)", r2_val)
cat(sprintf("   -> Model Fitted. R² = %.3f\n", r2_val))

# C. Prediction Grid
x_rng <- range(df_fit$PC1)
y_rng <- range(df_fit$PC2)
grid_data <- expand.grid(
  PC1 = seq(x_rng[1], x_rng[2], length.out = GRID_RES),
  PC2 = seq(y_rng[1], y_rng[2], length.out = GRID_RES)
)

# Predict
grid_data$Pred_Value <- predict(gam_model, newdata = grid_data)

# D. Masking (Convex Hull)
#    Remove predictions outside the actual data polygon
hull_idx <- chull(df_fit$PC1, df_fit$PC2)
hull_pts <- df_fit[hull_idx, c("PC1", "PC2")]
in_poly <- point.in.polygon(grid_data$PC1, grid_data$PC2, hull_pts$PC1, hull_pts$PC2)
grid_data$Pred_Value[in_poly == 0] <- NA

# E. Plotting
p_taxmpd <- ggplot() +
  # 1. Raster Layer (Smooth Surface)
  geom_raster(data = grid_data, aes(x = PC1, y = PC2, fill = Pred_Value), interpolate = TRUE) +
  
  # 2. Contour Lines
  geom_contour(data = grid_data, aes(x = PC1, y = PC2, z = Pred_Value), 
               colour = "grey50", linewidth = 0.3, bins = 6, na.rm = TRUE) +
  
  # 3. Color Scale
  #    Using linear gradient (Blue-White-Pink)
  scale_fill_gradientn(
    colors = CUSTOM_COLS,
    name = "SES Value",
    na.value = "transparent"
  ) +
  
  # 4. Hull Boundary (White)
  geom_polygon(data = hull_pts, aes(x = PC1, y = PC2), 
               fill = NA, colour = "white", linewidth = 1.2) +
  
  # 5. Theme & Layout
  labs(title = title_str, x = "PC1", y = "PC2") + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    
    # Force Square Aspect Ratio
    aspect.ratio = 1,
    
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title = element_text(color = "black", size = 12),
    axis.text  = element_text(color = "black", size = 10),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    
    # Legend Settings
    legend.position = "right",
    legend.key.height = unit(1.0, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.title = element_text(size = 10)
  )

# 5. Save Output
#-----------------------------------------------------------------------------
out_file <- "Fig4_d.pdf"

# Save as PDF (Square dimensions for single plot)
ggsave(out_file, p_taxmpd, width = 7, height = 7, dpi = 300)

cat(sprintf("\n>>> Success! Fig4_d saved to: %s\n", file.path(getwd(), out_file)))

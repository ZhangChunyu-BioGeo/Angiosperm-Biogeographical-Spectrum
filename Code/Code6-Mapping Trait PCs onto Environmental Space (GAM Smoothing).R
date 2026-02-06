#-----------------------------------------------------------------------------
# Code 6: Mapping Trait PCs onto Environmental Space (GAM Smoothing)
# 
# Description:
# 1. Performs PCA on Environmental variables (EnvPCA) to define the climate space.
# 2. Maps Trait PC1, PC2, and PC3 scores onto this space using GAM smoothing.
# 3. Visualizes the gradients using a contour/raster plot
#
# Inputs:
# - Global_Ecoregion_Environment.csv
# - PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv
#
# Outputs:
# - Fig2_e.pdf
#-----------------------------------------------------------------------------

# 0. Global Parameters
#-----------------------------------------------------------------------------
# [Visualization Settings]
SHOW_ARROWS <- FALSE       # Show environmental vectors?
GAM_K       <- 15          # GAM smoothing basis dimension
GRID_RES    <- 300         # Grid resolution (300x300)

# [Color Palette] (Green -> Cream -> Purple/Pink)
# Mapping: Min -> Neg_Mid -> Zero -> Pos_Mid -> Max
CUSTOM_COLS <- c("#427638", "#9DC86E", "#FDFBE9", "#D27EAC", "#8C1C54")

# 1. Setup and Data Loading
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("data.table", "ggplot2", "mgcv", "sp", "grid", "ggrepel", "scales", "patchwork")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(ggplot2)
library(mgcv)
library(sp)
library(grid)       
library(ggrepel)    
library(scales)     
library(patchwork)

cat(">>> 1. Loading and Preparing Data...\n")

# --- File Paths ---
env_path   <- "Global_Ecoregion_Environment.csv"
# Updated input filename as requested
trait_path <- "PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv"

df_env   <- fread(env_path)
df_trait <- fread(trait_path)

# --- Merge Data ---
# Select necessary columns. Note: Input file usually has PC1, PC2, PC3 headers.
trait_subset <- df_trait[, .(ECO_ID, Trait_PC1 = PC1, Trait_PC2 = PC2, Trait_PC3 = PC3)]
merged_data  <- merge(df_env, trait_subset, by = "ECO_ID")
merged_data  <- na.omit(merged_data)

cat(sprintf("   -> Analysis based on %d ecoregions.\n", nrow(merged_data)))

# 2. Environmental PCA (EnvPCA)
#-----------------------------------------------------------------------------
cat(">>> 2. Performing Environmental PCA...\n")

env_vars <- c("AMT", "MDTR", "TS", "AP", "PS", "AI", "PET", "SRAD", 
              "Wind", "AMTd", "APd", "Elev", "Slope", "Sand", "SOC", "PH")

pca_res <- prcomp(merged_data[, ..env_vars], scale. = TRUE)

# --- Axis Direction Adjustment (Consistency) ---
env_pc_names <- paste0("EnvPC", 1:ncol(pca_res$x))
colnames(pca_res$x) <- env_pc_names
colnames(pca_res$rotation) <- env_pc_names

# 1. Align EnvPC1 with AMT (Positive)
if (cor(pca_res$x[, "EnvPC1"], merged_data$AMT) < 0) {
  pca_res$x[, "EnvPC1"] <- -pca_res$x[, "EnvPC1"]
  pca_res$rotation[, "EnvPC1"] <- -pca_res$rotation[, "EnvPC1"]
}

# 2. Force Flip EnvPC2 (as per previous logic)
pca_res$x[, "EnvPC2"] <- -pca_res$x[, "EnvPC2"]
pca_res$rotation[, "EnvPC2"] <- -pca_res$rotation[, "EnvPC2"]

# --- Prepare Analysis Data ---
df_analysis <- cbind(merged_data, as.data.frame(pca_res$x[, 1:2]))

# Calculate Variance Explained
eig_vals <- pca_res$sdev^2
var_exp  <- round(eig_vals / sum(eig_vals) * 100, 1)
pc1_lab  <- paste0("EnvPC1 (", var_exp[1], "%)")
pc2_lab  <- paste0("EnvPC2 (", var_exp[2], "%)")

# Prepare Loadings (Arrows) - Only if enabled
if (SHOW_ARROWS) {
  df_loadings <- as.data.frame(pca_res$rotation[, 1:2])
  df_loadings$var <- rownames(df_loadings)
  
  # Auto-scaling logic
  ratio_x <- max(abs(df_analysis$EnvPC1)) / max(abs(df_loadings$EnvPC1))
  ratio_y <- max(abs(df_analysis$EnvPC2)) / max(abs(df_loadings$EnvPC2))
  scale_factor <- min(ratio_x, ratio_y) * 0.75 
  
  df_loadings$EnvPC1_scaled <- df_loadings$EnvPC1 * scale_factor
  df_loadings$EnvPC2_scaled <- df_loadings$EnvPC2 * scale_factor
}

# 3. Define Plotting Function (GAM Surface)
#-----------------------------------------------------------------------------
plot_gam_surface <- function(trait_col, trait_display_name) {
  
  cat(sprintf("   -> Processing: %s ... ", trait_display_name))
  
  # A. Prepare Data
  df_fit <- df_analysis[, c("EnvPC1", "EnvPC2", trait_col), with = FALSE]
  setnames(df_fit, trait_col, "Trait_Value")
  
  # B. GAM Fitting
  gam_model <- gam(Trait_Value ~ s(EnvPC1, EnvPC2, k = GAM_K), data = df_fit)
  
  # Extract Adj. R-squared
  r2_val <- summary(gam_model)$r.sq
  title_str <- sprintf("%s\n(Adj. R² = %.3f)", trait_display_name, r2_val)
  cat(sprintf("R2 = %.3f\n", r2_val))
  
  # C. Prediction Grid
  x_rng <- range(df_fit$EnvPC1)
  y_rng <- range(df_fit$EnvPC2)
  grid_data <- expand.grid(
    EnvPC1 = seq(x_rng[1], x_rng[2], length.out = GRID_RES),
    EnvPC2 = seq(y_rng[1], y_rng[2], length.out = GRID_RES)
  )
  grid_data$Pred_Value <- predict(gam_model, newdata = grid_data)
  
  # D. Convex Hull Clipping (Prevent Extrapolation)
  hull_idx <- chull(df_fit$EnvPC1, df_fit$EnvPC2)
  hull_pts <- df_fit[hull_idx, c("EnvPC1", "EnvPC2")]
  in_poly <- point.in.polygon(grid_data$EnvPC1, grid_data$EnvPC2, hull_pts$EnvPC1, hull_pts$EnvPC2)
  grid_data$Pred_Value[in_poly == 0] <- NA
  
  # E. Color Mapping Logic (5-Point Divergent)
  val_min <- min(grid_data$Pred_Value, na.rm=TRUE)
  val_max <- max(grid_data$Pred_Value, na.rm=TRUE)
  
  # Ensure #FDFBE9 aligns with 0
  rescale_vals <- scales::rescale(c(val_min, val_min/2, 0, val_max/2, val_max))
  
  # F. Visualization
  p <- ggplot() +
    # Layer 1: Raster Heatmap
    geom_raster(data = grid_data, aes(x = EnvPC1, y = EnvPC2, fill = Pred_Value), interpolate = TRUE) +
    
    # Layer 2: Contours
    geom_contour(data = grid_data, aes(x = EnvPC1, y = EnvPC2, z = Pred_Value), 
                 colour = "grey50", linewidth = 0.3, bins = 6, na.rm = TRUE) +
    
    # Layer 3: Color Scale
    scale_fill_gradientn(
      colors = CUSTOM_COLS,
      values = rescale_vals,
      name = "Value",
      na.value = "transparent"
    ) +
    
    # Layer 4: Hull Boundary
    geom_polygon(data = hull_pts, aes(x = EnvPC1, y = EnvPC2), 
                 fill = NA, colour = "white", linewidth = 1.2) +
    
    # Optional Arrows
    {if(SHOW_ARROWS) geom_segment(data = df_loadings, 
                                  aes(x = 0, y = 0, xend = EnvPC1_scaled, yend = EnvPC2_scaled),
                                  arrow = arrow(length = unit(0.2, "cm")), 
                                  colour = "black", linewidth = 0.4, alpha=0.7)} +
    {if(SHOW_ARROWS) geom_text_repel(data = df_loadings, 
                                     aes(x = EnvPC1_scaled, y = EnvPC2_scaled, label = var),
                                     colour = "black", size = 3, bg.color = "white", bg.r = 0.15)} +
    
    # Theme & Titles
    labs(title = title_str, x = pc1_lab, y = pc2_lab) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.title = element_text(color = "black"),
      axis.text  = element_text(color = "black"),
      axis.ticks = element_line(colour = "black"),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "right",
      legend.key.height = unit(1.0, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.title = element_blank()
    ) +
    coord_fixed()
  
  return(p)
}

# 4. Generate and Stitch Plots
#-----------------------------------------------------------------------------
cat(">>> 3. Generating Plots...\n")

p1 <- plot_gam_surface("Trait_PC1", "Trait PC1")
p2 <- plot_gam_surface("Trait_PC2", "Trait PC2")
p3 <- plot_gam_surface("Trait_PC3", "Trait PC3")

cat(">>> 4. Stitching and Saving...\n")

# Use patchwork for robust alignment
combined_plot <- p1 + p2 + p3 + 
  plot_layout(ncol = 3) 

print(combined_plot)

# Save as PDF
ggsave("Fig2_e.pdf", combined_plot, width = 18, height = 6)
cat(">>> Saved: Fig2_e.pdf\n")

#-----------------------------------------------------------------------------
# Code5: Data Process of Environmental PCA and Biplot Visualization (Fig 2d)
# 
# Description:
# 1. Data Preparation: Merges Environment data with Trait PCA results.
# 2. Bivariate Classification: Creates 'PinkGrn' classes based on Trait PC1 & PC2.
# 3. Size Sorting: Orders 'SizeCategory' factor based on SpeciesCount values.
# 4. Environmental PCA: Performs PCA on 16 environmental variables.
# 5. Visualization (Fig 2d): Generates a biplot with:
#    - Points colored by bivariate class.
#    - Points sized by Species Richness category.
#    - Environmental vectors (arrows).
#    - Custom Bivariate + Size legend layout.
#
# Inputs:
# - Global_Ecoregion_Environment.csv
# - PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv
#
# Outputs:
# - Fig2_d.pdf
# - Global_Ecoregion_Environment_EnvPCs.csv
#-----------------------------------------------------------------------------

# 0. Global Parameters (User Settings)
#-----------------------------------------------------------------------------
# [Control Parameters] 

# A. Axis Limits
xlim_manual <- NULL       
ylim_manual <- c(-5, 5.5) 

# B. Arrow Settings
arrow_scale_factor <- 1.2   # Scaling factor for arrow length
arrow_head_size    <- 0.35  # Arrow head length (cm)

# 1. Set Working Directory & Libraries
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("data.table", "factoextra", "ggplot2", "biscale", "cowplot", "grid")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(factoextra)
library(ggplot2)
library(biscale)
library(cowplot)
library(grid)

# 2. Data Loading & Preparation
#-----------------------------------------------------------------------------
# --- File Paths ---
env_data_path  <- "Global_Ecoregion_Environment.csv"
trait_pca_path <- "PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv"

# --- Load Data ---
cat("Loading data...\n")
df_env <- fread(env_data_path)
df_trait <- fread(trait_pca_path) 

# --- Merge ---
# Note: Including SpeciesCount (or SpeciesNumber) for sorting logic
# Ensure the column name matches your file (using SpeciesNumber based on Code 1 output)
if("SpeciesNumber" %in% names(df_trait) && !"SpeciesCount" %in% names(df_trait)) {
  setnames(df_trait, "SpeciesNumber", "SpeciesCount")
}

trait_subset <- df_trait[, .(ECO_ID, Trait_PC1 = PC1, Trait_PC2 = PC2, SizeCategory, SpeciesCount)]
merged_data <- merge(df_env, trait_subset, by = "ECO_ID")
merged_data <- na.omit(merged_data)

cat(sprintf("Analysis based on %d ecoregions.\n", nrow(merged_data)))

# --- Bivariate Classification (PinkGrn) ---
merged_data <- bi_class(merged_data, x = Trait_PC1, y = Trait_PC2, style = "quantile", dim = 3)
custom_biscale_pal <- bi_pal("PinkGrn", dim = 3, preview = FALSE)

# --- Size Category Sorting (Fixed Logic) ---
# Logic: Determine order by the minimum SpeciesCount in each category (Small -> Large)
cat("Sorting Size Categories using SpeciesCount numeric values...\n")

# Calculate representative value (min) for each category
size_rank <- merged_data[, .(MinCount = min(SpeciesCount)), by = SizeCategory]
# Sort ascending
ordered_levels <- size_rank[order(MinCount)]$SizeCategory

# Apply factor order
merged_data$SizeCategory <- factor(merged_data$SizeCategory, levels = ordered_levels)

cat("Category Order (Small -> Large):\n")
print(levels(merged_data$SizeCategory))

# 3. Perform Environmental PCA
#-----------------------------------------------------------------------------
env_vars <- c("AMT", "MDTR", "TS", "AP", "PS", "AI", "PET", "SRAD", 
              "Wind", "AMTd", "APd", "Elev", "Slope", "Sand", "SOC", "PH")

pca_res <- prcomp(merged_data[, ..env_vars], scale. = TRUE)

# --- Direction Adjustment ---
env_pc_names <- paste0("EnvPC", 1:ncol(pca_res$x))
colnames(pca_res$x) <- env_pc_names
colnames(pca_res$rotation) <- env_pc_names

# 1. Align EnvPC1 with AMT (Positive correlation)
if (cor(pca_res$x[, "EnvPC1"], merged_data$AMT) < 0) {
  cat("Flipping EnvPC1...\n")
  pca_res$x[, "EnvPC1"] <- -pca_res$x[, "EnvPC1"]
  pca_res$rotation[, "EnvPC1"] <- -pca_res$rotation[, "EnvPC1"]
}

# 2. Flip EnvPC2 (User Request)
cat("Flipping EnvPC2 direction...\n")
pca_res$x[, "EnvPC2"] <- -pca_res$x[, "EnvPC2"]
pca_res$rotation[, "EnvPC2"] <- -pca_res$rotation[, "EnvPC2"]

# 3. Scale Arrows
if (arrow_scale_factor != 1) {
  pca_res$rotation <- pca_res$rotation * arrow_scale_factor
}

# 4. Visualization Construction
#-----------------------------------------------------------------------------
eig_val <- get_eigenvalue(pca_res)
pc1_lab <- paste0("EnvPC1 (", round(eig_val$variance.percent[1], 1), "%)")
pc2_lab <- paste0("EnvPC2 (", round(eig_val$variance.percent[2], 1), "%)")

# Point sizes corresponding to Small -> Large categories
size_values <- c(2, 3.5, 5) 

# --- Step 1: Base Biplot ---
p_base <- fviz_pca_biplot(
  pca_res,
  axes = c(1, 2),
  
  # A. Points
  geom.ind = "point",
  pointshape = 21,             
  fill.ind = merged_data$bi_class, 
  col.ind = "white",           
  pointsize = merged_data$SizeCategory, 
  alpha.ind = 0.9,
  stroke = 0.4,                
  
  # B. Variables
  geom.var = c("arrow", "text"),
  col.var = "black",
  alpha.var = 0.8,
  repel = TRUE,                
  
  # C. Theme
  title = "",
  ggtheme = theme_minimal()
) +
  labs(x = pc1_lab, y = pc2_lab) +
  scale_fill_manual(values = custom_biscale_pal, name = "") +
  
  # D. Size Legend
  scale_size_manual(
    name = "Species Richness", 
    values = size_values,
    breaks = levels(merged_data$SizeCategory), # Ensure correct legend order
    guide = guide_legend(
      override.aes = list(
        shape = 21, 
        fill = "grey50",       
        color = "white",       
        stroke = 0.5
      )
    )
  ) +
  
  # E. Theme Optimization (Refined)
  theme(
    # Panel background & Border
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Axis Text
    axis.text = element_text(color = "black", size = 10),
    
    # Axis Titles (Plain)
    axis.title = element_text(size = 12), 
    
    # Ticks (0.2cm)
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm")
  )

# --- Step 1.5: Arrow Head Modification & Axis Limits ---
# Access ggplot layers to modify arrow head size
if (length(p_base$layers) >= 2) {
  if (inherits(p_base$layers[[2]]$geom, "GeomSegment")) {
    p_base$layers[[2]]$aes_params$arrow <- arrow(length = unit(arrow_head_size, "cm"))
  }
}

# Apply manual limits if set
if (!is.null(xlim_manual) || !is.null(ylim_manual)) {
  p_base <- p_base + coord_cartesian(xlim = xlim_manual, ylim = ylim_manual)
}

# --- Step 2: Legends Extraction & Assembly ---
# Create a dummy plot to extract the Size legend
p_legend_extract <- p_base + theme(
  legend.position = "right",
  legend.title = element_text(size = 10) 
)
size_legend <- get_legend(p_legend_extract)

# Generate Biscale Legend
bi_legend_custom <- bi_legend(
  pal = custom_biscale_pal,
  dim = 3,
  xlab = "Trait PC1",
  ylab = "Trait PC2",
  size = 9
) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 9) 
  )

# Stack legends vertically
legend_col <- plot_grid(
  bi_legend_custom,
  size_legend,
  ncol = 1,
  rel_heights = c(1, 0.8), 
  align = "v"
)

# Combine Main Plot + Legend Column
final_plot <- plot_grid(
  p_base + theme(legend.position = "none"), 
  legend_col,
  ncol = 2,
  rel_widths = c(1, 0.25) 
)

# 5. Output
#-----------------------------------------------------------------------------
print(final_plot)

# Save Plot
ggsave("Fig2_d.pdf", final_plot, width = 8, height = 6)
cat("  -> Saved: Fig2_d.pdf\n")

# Save EnvPC Results
output_data <- cbind(merged_data[, .(ECO_ID)], as.data.table(pca_res$x[, 1:3]))
fwrite(output_data, "Global_Ecoregion_Environment_EnvPCs.csv")
cat("  -> Saved: Global_Ecoregion_Environment_EnvPCs.csv\n")

cat("\nAnalysis Completed!\n")

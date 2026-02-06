#-----------------------------------------------------------------------------
# Code 23: Global Spatial Piecewise SEM (Fig 5)
#
# Description:
#   Performs Piecewise Structural Equation Modeling (pSEM) to disentangle 
#   the relationships between Environment, Phylogenetic History, Trait Composition, 
#   and Functional Divergence (TaxMPD_SES).
#
#   Uses 3D Cartesian coordinates (X, Y, Z) in the correlation 
#     structure (GLS) to account for global spatial autocorrelation.
#
# Inputs:
#   1. Diversity_Taxonomic_Traits_SES_Summary.csv (Response: TaxMPD_SES)
#   2. Global_Ecoregion_Environment_EnvPCs.csv (Predictors: EnvPC1, EnvPC2)
#   3. Global_Ecoregion_Phylogenetic.csv (Predictor: MPD - Phylogenetic History)
#   4. PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv (Mediators: PC1, PC2)
#   5. Global_Ecoregion_Coordinates.csv (Spatial: Lat, Lon)
#
# Outputs:
#   - Global_Spatial_SEM_Summary.txt (Model Fit & D-sep tests)
#   - Global_Spatial_SEM_Coefs.csv (Path Coefficients)
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

# Load Packages
pkgs <- c("piecewiseSEM", "nlme", "data.table", "stats", "dplyr")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
}

library(piecewiseSEM)
library(nlme)
library(data.table)
library(dplyr)

# 2. Load Data
#-----------------------------------------------------------------------------
cat("1. Loading datasets...\n")

# Define file paths
f_tax   <- "Diversity_Taxonomic_Traits_SES_Summary.csv"
f_env   <- "Global_Ecoregion_Environment_EnvPCs.csv"
f_phy   <- "Global_Ecoregion_Phylogenetic.csv"
f_trait <- "PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv"
f_coord <- "Global_Ecoregion_Coordinates.csv"

# Check existence
if(!all(file.exists(f_tax, f_env, f_phy, f_trait, f_coord))) {
  stop("Error: One or more input files are missing.")
}

# Read Data
df_tax_ses   <- fread(f_tax)
df_env       <- fread(f_env)
df_phylo     <- fread(f_phy)
df_trait_pca <- fread(f_trait)
df_coords    <- fread(f_coord)

# 3. Data Cleaning, Merging, and Standardization
#-----------------------------------------------------------------------------
cat("2. Processing Data...\n")

# Ensure ID columns are character type
df_tax_ses$ECO_ID   <- as.character(df_tax_ses$ECO_ID)
df_env$ECO_ID       <- as.character(df_env$ECO_ID)
df_phylo$ECO_ID     <- as.character(df_phylo$ECO_ID)
df_trait_pca$ECO_ID <- as.character(df_trait_pca$ECO_ID)
df_coords$ECO_ID    <- as.character(df_coords$ECO_ID)

# Merge all datasets by ECO_ID
data_list <- list(df_tax_ses, df_env, df_phylo, df_trait_pca, df_coords)
analysis_data <- Reduce(function(x, y) merge(x, y, by = "ECO_ID", all = FALSE), data_list)

# Remove NAs
analysis_data <- na.omit(analysis_data)

# [Standardization - CRITICAL FIX] 
# Scale continuous variables and force them to be numeric VECTORS.
# 'scale()' normally returns a matrix, which causes piecewiseSEM summary to crash.
cols_to_scale <- c("TaxMPD_SES", "EnvPC1", "EnvPC2", "PC1", "PC2", "MPD")

analysis_data_scaled <- analysis_data %>%
  mutate(across(all_of(cols_to_scale), ~ as.numeric(scale(.)))) %>% 
  as.data.frame() # piecewiseSEM requires a plain data.frame, not data.table

cat(sprintf("   - Dataset ready. N = %d ecoregions.\n", nrow(analysis_data_scaled)))

# 4. Convert Spherical Coordinates to 3D Cartesian (X, Y, Z)
#-----------------------------------------------------------------------------
cat("3. Converting Spherical Coordinates to 3D Cartesian (X, Y, Z)...\n")

R <- 6371 # Earth Radius (km)
rad <- pi / 180

lat_rad <- analysis_data_scaled$Lat * rad
lon_rad <- analysis_data_scaled$Lon * rad 

analysis_data_scaled$X <- R * cos(lat_rad) * cos(lon_rad)
analysis_data_scaled$Y <- R * cos(lat_rad) * sin(lon_rad)
analysis_data_scaled$Z <- R * sin(lat_rad)

cat("   - Conversion done.\n")

# 5. Define Spatial Correlation Structure
#-----------------------------------------------------------------------------
# Using Exponential correlation with nugget effect based on 3D coordinates
spatial_cor <- corExp(form = ~ X + Y + Z, nugget = TRUE)

# 6. Fit GLS Sub-models
#-----------------------------------------------------------------------------
cat("4. Fitting GLS sub-models with 3D spatial correlation...\n")
cat("   (This may take a minute due to spatial matrix inversion...)\n")

# [Model 1] Phylogenetic Structure (MPD)
# Hypothesis: Driven by Environment
model_mpd <- gls(
  MPD ~ EnvPC1 + EnvPC2,
  correlation = spatial_cor,
  data = analysis_data_scaled,
  method = "ML"
)

# [Model 2] Trait Composition PC1
# Hypothesis: Driven by Environment + Phylogeny
model_pc1 <- gls(
  PC1 ~ EnvPC1 + EnvPC2 + MPD,
  correlation = spatial_cor,
  data = analysis_data_scaled,
  method = "ML" 
)

# [Model 3] Trait Composition PC2
# Hypothesis: Driven by Environment + Phylogeny
model_pc2 <- gls(
  PC2 ~ EnvPC1 + EnvPC2 + MPD,
  correlation = spatial_cor,
  data = analysis_data_scaled,
  method = "ML"
)

# [Model 4] Functional Divergence (TaxMPD_SES)
# Hypothesis: The final response, driven by everything
model_tax <- gls(
  TaxMPD_SES ~ EnvPC1 + PC1 + PC2 + MPD,
  correlation = spatial_cor,
  data = analysis_data_scaled,
  method = "ML"
)

# 7. Assemble and Run Piecewise SEM
#-----------------------------------------------------------------------------
cat("5. Assembling and Testing Piecewise SEM...\n")

# Define the pSEM structure
sem_global <- psem(
  model_mpd,
  model_pc1,
  model_pc2,
  model_tax,
  PC1 %~~% PC2, 
  data = analysis_data_scaled
)

# Run Summary (Performs d-sep tests for missing paths)
summary_res <- tryCatch({
  summary(sem_global)
}, error = function(e) {
  cat("Error in summary(sem_global): ", e$message, "\n")
  return(NULL)
})

# 8. Save Results
#-----------------------------------------------------------------------------
cat("\n6. Saving Results...\n")

# A. Save Text Summary (Fit statistics & D-separation)
out_txt <- "Global_Spatial_SEM_Summary.txt"
if (!is.null(summary_res)) {
  sink(out_txt)
  print(summary_res)
  sink()
  cat("Summary saved to:", out_txt, "\n")
} else {
  cat("Skipped saving summary text due to calculation error.\n")
}

# B. Save Coefficients (Path strengths)
out_csv <- "Global_Spatial_SEM_Coefs.csv"

# Ensure we can extract coefficients even if summary failed on d-sep
coef_df <- tryCatch({
  coefs(sem_global)
}, error = function(e) {
  cat("Error extracting coefficients: ", e$message, "\n")
  return(NULL)
})

if (!is.null(coef_df)) {
  write.csv(coef_df, out_csv, row.names = FALSE)
  cat("Coefficients saved to:", out_csv, "\n")
  cat("\nPreview of Coefficients:\n")
  print(head(coef_df))
}

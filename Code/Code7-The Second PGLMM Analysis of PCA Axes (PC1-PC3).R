#-----------------------------------------------------------------------------
# Code 7: The Second PGLMM Analysis of PCA Axes (PC1-PC3) - Space Only
#
# Description:
# Integrated PGLMM Analysis Pipeline (Spatial Random Effect Only + Phylo Fixed)
#
# Model Structure:
#   Response ~ EnvPC1 + EnvPC2 + PD + MDT + MPD + Area + (1|Space)
#
# Key Changes:
# 1. Phylogenetic metrics (PD, MDT, MPD) are treated as Fixed Effects.
# 2. Phylogenetic Random Effect is REMOVED.
# 3. SpeciesCount is REMOVED; only Area is retained as a control.
# 4. REMOVED all Linear Model (LM) calculations and comparisons.
#
# Inputs:
#   1. PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv
#   2. Global_Ecoregion_Environment_EnvPCs.csv
#   3. Global_Ecoregion_Phylogenetic.csv
#   4. Global_Ecoregion_Coordinates.csv
#   5. Global_Ecoregion_Space_Similarity_Matrix_IDW.csv
#
# Outputs:
#   1. PGLMM_Coefficients_Standardized_PCs_Phylo_Metric.csv
#   2. PGLMM_Model_Statistics_PCs_Phylo_Metric.csv
#-----------------------------------------------------------------------------

# 0. Global Settings and Library Loading
#-----------------------------------------------------------------------------
options(Matrix.warnDeprecatedCoerce = 0)
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("data.table", "phyr", "Matrix", "stringr")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) install.packages(p)
}

library(data.table)
library(phyr)
library(Matrix)
library(stringr)

# 1. Load Data
#-----------------------------------------------------------------------------
cat("1. Reading Data...\n")

# Define Paths
trait_path  <- "PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv"
env_path    <- "Global_Ecoregion_Environment_EnvPCs.csv"
phylo_path  <- "Global_Ecoregion_Phylogenetic.csv"
coord_path  <- "Global_Ecoregion_Coordinates.csv"
space_path  <- "Global_Ecoregion_Space_Similarity_Matrix_IDW.csv"

# Read Data
df_traits    <- fread(trait_path)
df_env       <- fread(env_path, select = c("ECO_ID", "EnvPC1", "EnvPC2"))
df_phylo_var <- fread(phylo_path, select = c("ECO_ID", "PD", "MDT", "MPD"))
df_coords    <- fread(coord_path, select = c("ECO_ID", "Area"))
df_space_raw <- fread(space_path)

# 2. Data Alignment
#-----------------------------------------------------------------------------
cat("2. Aligning Data...\n")

# Standardize ID column type
df_traits$ECO_ID    <- as.character(df_traits$ECO_ID)
df_env$ECO_ID       <- as.character(df_env$ECO_ID)
df_phylo_var$ECO_ID <- as.character(df_phylo_var$ECO_ID)
df_coords$ECO_ID    <- as.character(df_coords$ECO_ID)
df_space_raw$ECO_ID <- as.character(df_space_raw[[1]]) 

# Find Common IDs
common_ids <- Reduce(intersect, list(
  df_traits$ECO_ID,
  df_env$ECO_ID,
  df_phylo_var$ECO_ID,
  df_coords$ECO_ID,
  df_space_raw$ECO_ID
))

cat(sprintf("   Common Ecoregions (N): %d\n", length(common_ids)))

if (length(common_ids) == 0) stop("Error: No common ECO_ID found!")

# Alignment Function
align_df <- function(df, ids) {
  df_sub <- df[match(ids, df$ECO_ID), ]
  return(na.omit(df_sub)) 
}

# Align datasets
df_traits_sub <- align_df(df_traits, common_ids)
final_ids     <- df_traits_sub$ECO_ID

# Re-align others based on final valid IDs
df_env_sub    <- df_env[match(final_ids, ECO_ID), ]
df_phylo_sub  <- df_phylo_var[match(final_ids, ECO_ID), ]
df_coords_sub <- df_coords[match(final_ids, ECO_ID), ]

# Spatial Matrix Slicing
cat("   Aligning spatial matrix...\n")
space_mat <- as.matrix(df_space_raw[, -1]) 
rownames(space_mat) <- as.character(df_space_raw[[1]])
colnames(space_mat) <- names(df_space_raw)[-1]
space_pd <- space_mat[final_ids, final_ids]

# 3. Variable Preparation and Standardization
#-----------------------------------------------------------------------------
cat("3. Preparing Variables...\n")

# Define Predictors
predictors_env   <- c("EnvPC1", "EnvPC2")
predictors_phylo <- c("PD", "MDT", "MPD")
predictors_ctrl  <- c("Area")

# Merge Fixed Effects and Scale
X_fixed_raw <- cbind(
  df_env_sub[, ..predictors_env],
  df_phylo_sub[, ..predictors_phylo],
  df_coords_sub[, ..predictors_ctrl]
)

# Standardization (Z-score)
X_fixed_scaled <- as.data.frame(scale(X_fixed_raw))

# Merge into Model Data
model_data <- cbind(df_traits_sub, X_fixed_scaled)
model_data$ECO_ID_Space <- as.factor(model_data$ECO_ID) # Random Effect Group

# Target Traits
trait_names <- c("PC1", "PC2", "PC3")

# 4. Analysis Loop: PGLMM Only (No LM)
#-----------------------------------------------------------------------------
results_coef_list <- list()
results_stat_list <- list()

cat("\n4. Starting PGLMM Analysis Loop (Space Random Only)...\n")

# Helper: Extract Variance Components (Space & Residual)
get_var_components_space_only <- function(mod) {
  # Residual Variance
  var_resid <- NA_real_
  if (!is.null(mod$s2resid)) {
    var_resid <- as.numeric(mod$s2resid)
  } else if (!is.null(mod$sigma)) {
    var_resid <- (mod$sigma)^2
  }
  
  # Space Variance
  var_space <- 0
  var_phylo <- 0 # Always 0 in this model configuration
  
  if (!is.null(mod$s2r)) {
    nm <- names(mod$s2r)
    idx_s <- grep("Space", nm)
    if (length(idx_s) > 0) var_space <- as.numeric(mod$s2r[idx_s[1]])
  } else if (!is.null(mod$ss)) {
    nm <- names(mod$ss)
    idx_s <- grep("Space", nm)
    if (length(idx_s) > 0) var_space <- as.numeric(mod$ss[idx_s[1]])^2
  }
  
  list(var_phylo = var_phylo, var_space = var_space, var_resid = var_resid)
}

# Helper: Get Variable Names
get_fixed_names <- function(mod, betas_length) {
  vars <- names(mod$B)
  if (is.null(vars) || all(vars == "")) {
    if (!is.null(colnames(mod$X))) vars <- colnames(mod$X)
    else vars <- paste0("X", seq_len(betas_length))
  }
  if (length(vars) != betas_length) vars <- paste0("X", seq_len(betas_length))
  vars
}

for (i in seq_along(trait_names)) {
  target_trait <- trait_names[i]
  cat(sprintf("[%d/%d] Fitting: %s ... ", i, length(trait_names), target_trait))
  
  # Predictor string: Env + Phylo + Area
  fixed_part <- paste(c(predictors_env, predictors_phylo, predictors_ctrl), collapse = " + ")
  
  # --- PGLMM (Space Random Effect Only) ---
  # Formula: Trait ~ Fixed + (1|Space)
  f_pglmm <- as.formula(paste(
    target_trait, "~", fixed_part, "+ (1|ECO_ID_Space)"
  ))
  
  tryCatch({
    mod <- pglmm(
      f_pglmm,
      data      = model_data,
      family    = "gaussian",
      cov_ranef = list(ECO_ID_Space = space_pd), # Space Matrix Only
      REML      = TRUE,
      verbose   = FALSE
    )
    
    # 1. Extract Fixed Coefficients
    if (!is.null(mod$B)) {
      betas <- as.numeric(mod$B)
      ses   <- as.numeric(mod$B.se)
      pvals <- as.numeric(mod$B.pvalue)
      vars  <- get_fixed_names(mod, length(betas))
      
      coef_dt <- data.table(
        Trait     = target_trait,
        Variable  = vars,
        Std_Beta  = betas,
        Std_Error = ses,
        P_Value   = pvals
      )
      results_coef_list[[i]] <- coef_dt
    } else {
      # Fallback
      betas <- as.numeric(coef(mod))
      results_coef_list[[i]] <- data.table(Trait = target_trait, Variable = names(coef(mod)), Std_Beta = betas, Std_Error = NA, P_Value = NA)
    }
    
    # 2. Extract Variance Components
    vc <- get_var_components_space_only(mod)
    var_space <- vc$var_space
    var_resid <- vc$var_resid
    var_phylo <- 0 
    
    # 3. Calculate R2 (Nakagawa)
    r2_m <- NA_real_; r2_c <- NA_real_; var_fixed <- NA_real_
    
    if (!is.null(mod$X) && !is.null(mod$B)) {
      fixed_pred <- as.vector(mod$X %*% mod$B)
      var_fixed  <- var(fixed_pred)
      var_total  <- var_fixed + var_space + var_resid 
      
      if (is.finite(var_total) && var_total > 0) {
        r2_m <- var_fixed / var_total
        r2_c <- (var_fixed + var_space) / var_total # Only adding Space variance
      }
    }
    
    # 4. Calculate ICC (Space Only)
    icc_space <- NA_real_; icc_phylo <- 0
    if (is.finite(var_resid)) {
      total_random_resid <- var_space + var_resid
      if (total_random_resid > 0) {
        icc_space <- var_space / total_random_resid
      } else {
        icc_space <- 0
      }
    }
    
    # 5. Store Statistics
    results_stat_list[[i]] <- data.table(
      Trait     = target_trait,
      LogLik    = as.numeric(mod$logLik),
      Var_Fixed = var_fixed,
      Var_Phylo = 0,             
      Var_Space = var_space,
      Var_Resid = var_resid,
      ICC_Phylo = 0,             
      ICC_Space = icc_space,
      R2_m      = r2_m,
      R2_c      = r2_c
    )
    
    cat(sprintf("Done. R2m=%.3f | R2c=%.3f\n", r2_m, r2_c))
    
  }, error = function(e) {
    cat(sprintf("FAILED! Msg: %s\n", e$message))
    results_stat_list[[i]] <- data.table(Trait=target_trait, LogLik=NA, Var_Fixed=NA, Var_Phylo=NA, Var_Space=NA, Var_Resid=NA, ICC_Phylo=NA, ICC_Space=NA, R2_m=NA, R2_c=NA)
  })
}

# 5. Export Results
#-----------------------------------------------------------------------------
cat("\n5. Exporting Results...\n")

# (1) Coefficient Table
final_coefs <- rbindlist(results_coef_list, fill = TRUE)
if (nrow(final_coefs) > 0) {
  final_coefs[, Significance := "ns"]
  final_coefs[!is.na(P_Value) & P_Value < 0.05,  Significance := "*"]
  final_coefs[!is.na(P_Value) & P_Value < 0.01,  Significance := "**"]
  final_coefs[!is.na(P_Value) & P_Value < 0.001, Significance := "***"]
  
  setcolorder(final_coefs, c("Trait", "Variable", "Std_Beta", "Std_Error", "P_Value", "Significance"))
  
  outfile1 <- "PGLMM_Coefficients_Standardized_PCs_Phylo_Metric.csv"
  fwrite(final_coefs, outfile1)
  cat(sprintf("  -> Saved: %s\n", outfile1))
}

# (2) Model Statistics Table
final_stats <- rbindlist(results_stat_list, fill = TRUE)
cols_stats <- c("Trait", "R2_m", "R2_c", "ICC_Phylo", "ICC_Space",
                "Var_Fixed", "Var_Phylo", "Var_Space", "Var_Resid", "LogLik")
final_stats <- final_stats[, ..cols_stats]

outfile2 <- "PGLMM_Model_Statistics_PCs_Phylo_Metric.csv"
fwrite(final_stats, outfile2)
cat(sprintf("  -> Saved: %s\n", outfile2))

cat("\nAnalysis Completed.\n")

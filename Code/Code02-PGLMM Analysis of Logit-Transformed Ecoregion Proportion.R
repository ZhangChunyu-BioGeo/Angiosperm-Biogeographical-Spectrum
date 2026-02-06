#-----------------------------------------------------------------------------
# Code 2: PGLMM Analysis of Logit-Transformed Ecoregion Proportions
#
# Description:
# 1. Data Preparation: Aligns trait, environment, phylogenetic, and spatial data.
#    - Input traits are already Logit-transformed (from Code 1).
# 2. PGLMM Fitting: Fits Phylogenetic Generalized Linear Mixed Models for each trait
#    using the 'phyr' package.
#    - Fixed Effects: Standardized environmental variables + Control variables (Area, SpeciesNumber).
#    - Random Effects: Phylogenetic signal and Spatial signal.
# 3. Result Extraction: Extracts standardized coefficients, variance components, 
#    and calculates R2_marginal (fixed effects) and R2_conditional (fixed + random).
#
# Output Files:
# 1. PGLMM_Coefficients_Standardized_Logit_Proportion.csv
# 2. PGLMM_Model_Statistics_Logit_Proportion.csv
#
# Last Update: 2026/02/06
#-----------------------------------------------------------------------------

# 0. Global Settings and Library Loading
#-----------------------------------------------------------------------------
options(Matrix.warnDeprecatedCoerce = 0)
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("data.table", "phyr", "Matrix", "stringr")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(phyr)
library(Matrix)
library(stringr)

# 1. Load Data
#-----------------------------------------------------------------------------
cat("1. Reading Data...\n")
# Input now points directly to the Logit-transformed file
trait_path  <- "Trait_Logit_Proportion_Ecoregion_No_Imputation_30Greater.csv"
env_path    <- "Global_Ecoregion_Environment.csv"
coord_path  <- "Global_Ecoregion_Coordinates.csv"
phylo_path  <- "Global_Ecoregion_Phylogenetic_Similarity_Matrix_PhyloSor.csv"
space_path  <- "Global_Ecoregion_Space_Similarity_Matrix_IDW.csv"

df_traits    <- fread(trait_path)
df_env       <- fread(env_path)
df_coords    <- fread(coord_path)
df_phylo_raw <- fread(phylo_path)
df_space_raw <- fread(space_path)

# 2. Data Alignment
#-----------------------------------------------------------------------------
cat("2. Aligning Data...\n")
ids_trait <- as.character(df_traits$ECO_ID)
ids_env   <- as.character(df_env$ECO_ID)
ids_coord <- as.character(df_coords$ECO_ID)
ids_phylo <- as.character(df_phylo_raw[[1]])
ids_space <- as.character(df_space_raw[[1]])

common_ids <- Reduce(intersect, list(ids_trait, ids_env, ids_coord, ids_phylo, ids_space))
cat(sprintf("   Common Ecoregions (N): %d\n", length(common_ids)))

if(length(common_ids) == 0) stop("Error: No common ECO_ID found!")

# Helper function to align and remove NAs
align_df <- function(df, ids) {
  df_sub <- df[match(ids, as.character(ECO_ID)), ]
  return(na.omit(df_sub))
}

df_traits_sub <- align_df(df_traits, common_ids)
final_ids     <- as.character(df_traits_sub$ECO_ID)

df_env_sub    <- df_env[match(final_ids,   as.character(ECO_ID)), ]
df_coords_sub <- df_coords[match(final_ids, as.character(ECO_ID)), ]

# Matrix slicing and alignment
cat("   Aligning matrices...\n")
phylo_mat <- as.matrix(df_phylo_raw[, -1])
rownames(phylo_mat) <- as.character(df_phylo_raw[[1]])
colnames(phylo_mat) <- names(df_phylo_raw)[-1]
phylo_pd <- phylo_mat[final_ids, final_ids]

space_mat <- as.matrix(df_space_raw[, -1])
rownames(space_mat) <- as.character(df_space_raw[[1]])
colnames(space_mat) <- names(df_space_raw)[-1]
space_pd <- space_mat[final_ids, final_ids]

# 3. Variable Preparation
#-----------------------------------------------------------------------------
cat("3. Preparing Variables (Standardization)...\n")

# (1) Environmental Variables
env_vars <- c("AMT", "MDTR", "AP", "PS", "SRAD", 
              "Wind", "AMTd", "APd", "Elev", "Slope", "Sand", "SOC", "PH")

# (2) Control Variables
# Note: Using SpeciesNumber from the input file
control_vars_raw <- data.frame(
  SpeciesNumber = df_traits_sub$SpeciesNumber,
  Area          = df_coords_sub$Area
)

# (3) Construct and standardize fixed effects matrix
X_fixed_raw    <- cbind(df_env_sub[, ..env_vars], control_vars_raw)
X_fixed_scaled <- as.data.frame(scale(X_fixed_raw))

# (4) Merge data for modeling (Traits are already Logit transformed)
model_data <- cbind(df_traits_sub, X_fixed_scaled)
model_data$ECO_ID_Phylo <- as.factor(model_data$ECO_ID)
model_data$ECO_ID_Space <- as.factor(model_data$ECO_ID)

# Target Trait Names
trait_names <- c("GF_H", "LC_A", "PH_D", "RS_T", "PT_U", "LS_C", "LM_L", "LV_P", 
                 "SS_B", "IF_I", "OP_H", "FS_Z", "PF_C", "FT_D", "FD_D")

# 4. Analysis Loop: PGLMM Fitting
#-----------------------------------------------------------------------------
results_coef_list <- list()
results_stat_list <- list()

cat("\n4. Starting PGLMM Analysis Loop...\n")

# Helper: Extract variance components from pglmm object
get_var_components_pglmm <- function(mod) {
  # Residual variance
  var_resid <- NA_real_
  if(!is.null(mod$s2resid)) {
    var_resid <- as.numeric(mod$s2resid)
  } else if (!is.null(mod$sigma)) {
    var_resid <- (mod$sigma)^2
  }
  
  # Random effect variances
  var_phylo <- 0
  var_space <- 0
  
  if (!is.null(mod$s2r)) {
    nm <- names(mod$s2r)
    if (!is.null(nm)) {
      idx_p <- grep("Phylo", nm)
      if (length(idx_p) > 0) var_phylo <- as.numeric(mod$s2r[idx_p[1]])
      idx_s <- grep("Space", nm)
      if (length(idx_s) > 0) var_space <- as.numeric(mod$s2r[idx_s[1]])
    }
  } else if (!is.null(mod$ss)) {
    nm <- names(mod$ss)
    if (!is.null(nm)) {
      idx_p <- grep("Phylo", nm)
      if (length(idx_p) > 0) var_phylo <- as.numeric(mod$ss[idx_p[1]])^2
      idx_s <- grep("Space", nm)
      if (length(idx_s) > 0) var_space <- as.numeric(mod$ss[idx_s[1]])^2
    }
  }
  
  list(var_phylo = var_phylo, var_space = var_space, var_resid = var_resid)
}

# Helper: Safely retrieve fixed effect names
get_fixed_names <- function(mod, betas_length) {
  vars <- names(mod$B)
  if (is.null(vars) || all(vars == "")) {
    if (!is.null(colnames(mod$X))) {
      vars <- colnames(mod$X)
    } else {
      vars <- paste0("X", seq_len(betas_length))
    }
  }
  if (length(vars) != betas_length) {
    vars <- paste0("X", seq_len(betas_length))
  }
  vars
}

for (i in seq_along(trait_names)) {
  target_trait <- trait_names[i]
  cat(sprintf("[%d/%d] Fitting: %s ... ", i, length(trait_names), target_trait))
  
  # Define predictors (Environment + Controls)
  predictors <- paste(c(env_vars, "SpeciesNumber", "Area"), collapse = " + ")
  
  # Define Formula: Logit(Trait) ~ Env + Controls + (1|Phylo) + (1|Space)
  f_full <- as.formula(paste(
    target_trait, "~", predictors,
    "+ (1|ECO_ID_Phylo) + (1|ECO_ID_Space)"
  ))
  
  tryCatch({
    mod <- pglmm(
      f_full,
      data      = model_data,
      family    = "gaussian",
      cov_ranef = list(ECO_ID_Phylo = phylo_pd, ECO_ID_Space = space_pd),
      REML      = TRUE,
      verbose   = FALSE
    )
    
    # 1. Extract Fixed Effect Coefficients
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
      # Fallback method
      betas <- as.numeric(coef(mod))
      vars  <- names(coef(mod))
      if (length(vars) != length(betas)) vars <- paste0("X", seq_along(betas))
      
      results_coef_list[[i]] <- data.table(
        Trait     = target_trait,
        Variable  = vars,
        Std_Beta  = betas,
        Std_Error = NA_real_,
        P_Value   = NA_real_
      )
    }
    
    # 2. Extract Variance Components
    vc <- get_var_components_pglmm(mod)
    var_phylo <- vc$var_phylo
    var_space <- vc$var_space
    var_resid <- vc$var_resid
    
    # 3. Calculate R2m and R2c (Nakagawa et al.)
    r2_m      <- NA_real_
    r2_c      <- NA_real_
    var_fixed <- NA_real_
    
    if (!is.null(mod$X) && !is.null(mod$B)) {
      fixed_pred <- as.vector(mod$X %*% mod$B)
      var_fixed  <- var(fixed_pred)
      
      var_total <- var_fixed + var_phylo + var_space + var_resid
      
      if (is.finite(var_total) && var_total > 0) {
        r2_m <- var_fixed / var_total
        r2_c <- (var_fixed + var_phylo + var_space) / var_total
      }
    }
    
    # 4. Calculate ICC
    icc_phylo <- NA_real_
    icc_space <- NA_real_
    if (is.finite(var_resid)) {
      total_random_resid <- var_phylo + var_space + var_resid
      if (total_random_resid > 0) {
        icc_phylo <- var_phylo / total_random_resid
        icc_space <- var_space / total_random_resid
      } else {
        icc_phylo <- 0
        icc_space <- 0
      }
    }
    
    # 5. Store Statistics
    results_stat_list[[i]] <- data.table(
      Trait     = target_trait,
      LogLik    = as.numeric(mod$logLik),
      Var_Fixed = var_fixed,
      Var_Phylo = var_phylo,
      Var_Space = var_space,
      Var_Resid = var_resid,
      ICC_Phylo = icc_phylo,
      ICC_Space = icc_space,
      R2_m      = r2_m,
      R2_c      = r2_c
    )
    
    cat(sprintf("Done. R2m=%.3f | R2c=%.3f\n", r2_m, r2_c))
    
  }, error = function(e) {
    cat(sprintf("FAILED! Msg: %s\n", e$message))
    results_stat_list[[i]] <- data.table(
      Trait = target_trait, LogLik = NA, Var_Fixed = NA, Var_Phylo = NA, 
      Var_Space = NA, Var_Resid = NA, ICC_Phylo = NA, ICC_Space = NA, 
      R2_m = NA, R2_c = NA
    )
  })
}

# 5. Export Results
#-----------------------------------------------------------------------------
cat("\n5. Exporting Results...\n")

# (1) Coefficient Table
final_coefs <- rbindlist(results_coef_list, fill = TRUE)

if (nrow(final_coefs) > 0) {
  # Add significance indicators
  final_coefs[, Significance := "ns"]
  final_coefs[!is.na(P_Value) & P_Value < 0.05,  Significance := "*"]
  final_coefs[!is.na(P_Value) & P_Value < 0.01,  Significance := "**"]
  final_coefs[!is.na(P_Value) & P_Value < 0.001, Significance := "***"]
  
  setcolorder(final_coefs, c("Trait", "Variable", "Std_Beta", "Std_Error", "P_Value", "Significance"))
  
  fwrite(final_coefs, "PGLMM_Coefficients_Standardized_Logit_Proportion.csv")
  cat("  -> Saved: PGLMM_Coefficients_Standardized_Logit_Proportion.csv\n")
}

# (2) Model Statistics Table
final_stats <- rbindlist(results_stat_list, fill = TRUE)
cols_stats <- c("Trait", "R2_m", "R2_c", "ICC_Phylo", "ICC_Space", 
                "Var_Fixed", "Var_Phylo", "Var_Space", "Var_Resid", "LogLik")
final_stats <- final_stats[, ..cols_stats]
fwrite(final_stats, "PGLMM_Model_Statistics_Logit_Proportion.csv")
cat("  -> Saved: PGLMM_Model_Statistics_Logit_Proportion.csv\n")

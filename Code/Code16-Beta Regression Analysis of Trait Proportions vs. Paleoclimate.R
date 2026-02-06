#-----------------------------------------------------------------------------
# Code 16: Beta Regression Analysis of Trait Proportions vs. Paleoclimate
#
# Description:
#   Uses Beta Regression to test the relationship between the proportions of 
#   30 trait states and Global Mean Surface Temperature (GMST).
#   Beta regression is suitable for dependent variables bounded in (0, 1).
#
# IMPORTANT NOTE:
#   The "Paleoclimate_Data_1Ma_Bins.csv" file is NOT included in this repository. 
#   Please obtain the PhanDA GMST data from Emily J. Judd et al.:
#   Source: https://github.com/EJJudd/PhanDA
#
# Inputs:
#   1. Evolution_SCM_Summary.csv (From Code 15)
#   2. Paleoclimate_Data_1Ma_Bins.csv (External Data)
#
# Outputs:
#   - GMST_Proportion_Coefficients_1MaBins.csv
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

if (!requireNamespace("betareg", quietly = TRUE)) install.packages("betareg")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

library(betareg)
library(dplyr)
library(tidyr)
library(readr)

# 2. Parameters and Paths
#-----------------------------------------------------------------------------
scm_file <- "Evolution_SCM_Summary.csv"
climate_file <- "Paleoclimate_Data_1Ma_Bins.csv"
output_file <- "GMST_Proportion_Coefficients_1MaBins.csv"

# Time window: Analyze only the last 100 Ma
MAX_AGE_BIN <- 100

# 3. Data Loading and Preparation
#-----------------------------------------------------------------------------
cat("1. Loading and Merging Data...\n")

if (!file.exists(climate_file)) {
  stop(paste("Error: Paleoclimate data not found.",
             "Please download GMST data from https://github.com/EJJudd/PhanDA",
             "and save as 'Paleoclimate_Data_1Ma_Bins.csv'."))
}

# Read SCM Summary
df_scm <- read_csv(scm_file, show_col_types = FALSE)

# Read Paleoclimate Data
df_climate <- read_csv(climate_file, show_col_types = FALSE) %>%
  dplyr::select(bin_id, GMST_50)

# Merge datasets
df_full <- inner_join(df_scm, df_climate, by = "bin_id")

# 4. Data Filtering and Transformation
#-----------------------------------------------------------------------------
cat("2. Filtering and Transforming Data (0-100 Ma)...\n")

# Filter time window
df_filtered <- df_full %>%
  filter(bin_id <= MAX_AGE_BIN)

# Pivot to Long Format for looping
df_long <- df_filtered %>%
  pivot_longer(
    cols = c("prop_State1_mean", "prop_State2_mean"),
    names_to = "Proportion_Type",
    values_to = "Proportion"
  ) %>%
  # Create unique label for Trait-State combination
  mutate(Trait_Type_Label = ifelse(
    Proportion_Type == "prop_State1_mean",
    State1_Label,
    State2_Label
  )) %>%
  # Transformation for Beta Regression:
  # Proportions must be strictly in (0, 1). 
  # Apply Smithson & Verkuilen (2006) transformation: y' = (y*(n-1) + 0.5)/n
  mutate(Proportion_Adj = (Proportion * (n_sims - 1) + 0.5) / n_sims) %>%
  dplyr::select(Trait, Trait_Type_Label, bin_id, Proportion_Adj, GMST_50)

# Get list of all trait-state combinations
unique_trait_types <- unique(df_long$Trait_Type_Label)
cat(sprintf("   - Starting Beta Regression for %d trait-state combinations.\n", length(unique_trait_types)))

# 5. Beta Regression Loop
#-----------------------------------------------------------------------------
cat("3. Executing Beta Regression...\n")

results_list <- list()
pb <- txtProgressBar(min = 0, max = length(unique_trait_types), style = 3)
counter <- 0

for (trait_type in unique_trait_types) {
  
  # Subset data
  df_subset <- df_long %>%
    filter(Trait_Type_Label == trait_type)
  
  # Fit Model
  fit <- tryCatch({
    betareg(Proportion_Adj ~ GMST_50, data = df_subset)
  }, error = function(e) {
    cat(sprintf("\nWarning: Fit failed for %s - %s\n", trait_type, e$message))
    return(NULL)
  })
  
  # Extract Coefficients if successful
  if (!is.null(fit)) {
    summary_fit <- summary(fit)
    
    # Extract info for GMST_50
    # summary_fit$coefficients$mean is a matrix, 2nd row is usually the predictor
    coef_info <- summary_fit$coefficients$mean["GMST_50", ]
    
    results_list[[trait_type]] <- data.frame(
      Trait_Type = trait_type,
      Coefficient = coef_info["Estimate"],
      Std_Error = coef_info["Std. Error"],
      Z_Value = coef_info["z value"],
      P_Value = coef_info["Pr(>|z|)"],
      Significance = dplyr::case_when(
        coef_info["Pr(>|z|)"] < 0.001 ~ "***",
        coef_info["Pr(>|z|)"] < 0.01  ~ "**",
        coef_info["Pr(>|z|)"] < 0.05  ~ "*",
        TRUE                         ~ "ns"
      ),
      stringsAsFactors = FALSE
    )
  }
  
  counter <- counter + 1
  setTxtProgressBar(pb, counter)
}
close(pb)

# 6. Aggregating and Saving Results
#-----------------------------------------------------------------------------
cat("\n4. Saving Results...\n")

if (length(results_list) > 0) {
  final_results_df <- bind_rows(results_list) %>%
    # Sort by P-value (most significant first)
    arrange(P_Value)
  
  write_csv(final_results_df, output_file)
  
  cat("Analysis Completed! Results saved to:\n", output_file, "\n\n")
  print(as.data.frame(final_results_df))
  
} else {
  cat("\nError: No valid regression results generated.\n")
}

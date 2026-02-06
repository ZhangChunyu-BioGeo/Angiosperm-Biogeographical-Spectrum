#-----------------------------------------------------------------------------
# Code 17: Data Integration of Trait Evolution Metrics
#
# Description:
#   Aggregates various evolutionary metrics into a single summary table for 
#   all 30 trait states (15 traits x 2 states).
#
# Inputs:
#   1. GMST_Proportion_Coefficients_1MaBins.csv (Beta Regression Coefs)
#   2. Phylogenetic_Signal_Summary.csv (Fritz's D)
#   3. Phylogenetic_consenTRAIT.csv (Evolutionary Depth)
#   4. Evolution_Models_Summary.csv (Transition Rates)
#   5. Loading_Trait_Logit_Proportion_Ecoregion_PCA_No_Imputation.csv (PCA Loadings)
#
# Outputs:
#   - Trait_Evolution_Summary.csv
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")
cat("Current Working Directory:", getwd(), "\n")

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(dplyr)
library(readr)
library(tidyr)

# 2. Load and Process Beta Regression Coefficients (Skeleton)
#-----------------------------------------------------------------------------
cat("1. Loading Paleoclimate Response Data (GMST Coefficients)...\n")
# This file contains 30 unique Trait_Type rows
df_gmst <- read_csv("GMST_Proportion_Coefficients_1MaBins.csv", show_col_types = FALSE) %>%
  select(
    Trait_Type, 
    GMST_Coef = Coefficient, 
    GMST_Coef_sd = Std_Error, 
    GMST_Coef_p = P_Value
  )

# Extract Base Trait (first 2 characters) for subsequent joining
df_base <- df_gmst %>%
  mutate(Trait = substr(Trait_Type, 1, 2))

# 3. Load and Merge Fritz's D (Phylogenetic Signal)
#-----------------------------------------------------------------------------
cat("2. Loading Fritz's D Data...\n")
# D values are calculated at the Trait level
df_fritz <- read_csv("Phylogenetic_Signal_Summary.csv", show_col_types = FALSE) %>%
  select(Trait, Fritz_D = Mean_D_Value, Fritz_D_sd = SD_D_Value, State_1 = State_0_Ref, State_2 = State_1_Alt)

# Join by Trait (Both states inherit the same D value)
df_base <- df_base %>%
  left_join(df_fritz %>% select(Trait, Fritz_D, Fritz_D_sd), by = "Trait")

# 4. Load and Merge consenTRAIT (Evolutionary Depth)
#-----------------------------------------------------------------------------
cat("3. Loading consenTRAIT Depth Data...\n")
# Modification:
# We assume Mean_Depth represents the divergence time depth for the trait.
# We select only Trait to join one-to-many, ignoring Target_State.
df_depth <- read_csv("Phylogenetic_consenTRAIT.csv", show_col_types = FALSE) %>%
  select(Trait, consenTRAIT_Depth = Mean_Depth, consenTRAIT_p = P_Value)

# Join by Trait
# Result: Both states (e.g., GF_W and GF_H) receive the depth value for GF
df_base <- df_base %>%
  left_join(df_depth, by = "Trait")

# 5. Process corHMM Transition Rates (Rate 1->2 and 2->1)
#-----------------------------------------------------------------------------
cat("4. Loading and Mapping corHMM Transition Rates...\n")

# Read model results
df_models <- read_csv("Evolution_Models_Summary.csv", show_col_types = FALSE)

# Read state definitions (to identify State 1 and State 2 names)
# Note: Ensure columns match those in Phylogenetic_Signal_Summary.csv or hardcode if needed
df_states_ref <- df_fritz %>% select(Trait, State_1, State_2)

# Combine model results with state definitions
df_rates_prep <- df_models %>%
  inner_join(df_states_ref, by = "Trait")

# Logic Split:
# Part A: If Trait_Type is State 1, we want the rate of leaving it (Rate_1_to_2)
df_rate_s1 <- df_rates_prep %>%
  select(Trait_Type = State_1, Transition_Rate = Rate_1_to_2)

# Part B: If Trait_Type is State 2, we want the rate of leaving it (Rate_2_to_1)
df_rate_s2 <- df_rates_prep %>%
  select(Trait_Type = State_2, Transition_Rate = Rate_2_to_1)

# Bind rows
df_rates_final <- bind_rows(df_rate_s1, df_rate_s2)

# Join to main table
df_base <- df_base %>%
  left_join(df_rates_final, by = "Trait_Type")

# 6. Load and Merge PCA Loadings (PC1, PC2, PC3)
#-----------------------------------------------------------------------------
cat("5. Loading PCA Loading Data...\n")
# The PCA file 'Trait' column corresponds to 'Trait_Type' here
df_pca <- read_csv("Loading_Trait_Logit_Proportion_Ecoregion_PCA_No_Imputation.csv", show_col_types = FALSE) %>%
  select(Trait_Type = Trait, PC1_Loading = PC1, PC2_Loading = PC2, PC3_Loading = PC3)

# Join to main table
df_base <- df_base %>%
  left_join(df_pca, by = "Trait_Type")

# 7. Final Organization and Saving
#-----------------------------------------------------------------------------
cat("6. Finalizing Table...\n")

final_df <- df_base %>%
  select(
    Trait,              # Base Trait Name (e.g., GF)
    Trait_Type,         # Trait State Name (e.g., GF_H)
    
    # 1. Paleoclimate Response
    GMST_Coef,
    GMST_Coef_sd,
    GMST_Coef_p,
    
    # 2. Phylogenetic Signal (Trait Level)
    Fritz_D,
    Fritz_D_sd,
    
    # 3. Evolutionary Depth (Trait Level)
    consenTRAIT_Depth,
    consenTRAIT_p,
    
    # 4. Transition Rate (Type Level: Rate LEAVING this state)
    Transition_Rate,
    
    # 5. Spatial Patterns (Type Level)
    PC1_Loading,
    PC2_Loading,
    PC3_Loading
  ) %>%
  arrange(Trait, Trait_Type)

# Check row count
cat("   - Final Row Count:", nrow(final_df), "\n")

# Save
output_path <- "Trait_Evolution_Summary.csv"
write_csv(final_df, output_path)

cat("\n==============================================================\n")
cat("Integration Completed! File saved to:\n", output_path, "\n")
cat("==============================================================\n")

# Preview
print(head(final_df))

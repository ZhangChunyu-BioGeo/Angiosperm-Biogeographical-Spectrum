#-----------------------------------------------------------------------------
# Code 12: Evolutionary Model Fitting (ER vs ARD)
#
# Description:
# Fits two standard models of binary character evolution for 15 functional traits
# using the `corHMM` package:
#   1. ER (Equal Rates): Transition rates between states are equal (q01 = q10).
#   2. ARD (All Rates Different): Transition rates are asymmetric (q01 != q10).
#
# Inputs:
#   1. Tree: GBOTB_Flora_GlobalID.tre
#   2. Data: Trait_Data_No_Imputation_QC_Ref.csv
#
# Outputs:
#   - ./Flora_Evolution_Models/{Trait}/ER_fit.rds
#   - ./Flora_Evolution_Models/{Trait}/ARD_fit.rds
#   - ./Flora_Evolution_Models/{Trait}/Model_Comparison_{Trait}.csv
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("corHMM", quietly = TRUE)) install.packages("corHMM")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ape)
library(corHMM)
library(dplyr)

# 2. Parameter Settings
#-----------------------------------------------------------------------------

# [Control Parameter: Select Trait to Analyze]
# Options: "GF", "LC", "PH", "RS", "PT", "LS", "LM", "LV", 
#          "SS", "IF", "OP", "FS", "PF", "FT", "FD"
TRAIT_TO_ANALYZE <- "GF" 

# 3. Trait Definition Dictionary (Defines State Order: 1 -> 2)
#-----------------------------------------------------------------------------
TRAIT_DEFINITIONS <- list(
  "GF" = c("GF_W", "GF_H"),
  "LC" = c("LC_P", "LC_A"),
  "PH" = c("PH_E", "PH_D"),
  "RS" = c("RS_T", "RS_F"),
  "PT" = c("PT_U", "PT_M"),
  "LS" = c("LS_S", "LS_C"),
  "LM" = c("LM_E", "LM_L"),
  "LV" = c("LV_R", "LV_P"),
  "SS" = c("SS_B", "SS_U"),
  "IF" = c("IF_D", "IF_I"),
  "OP" = c("OP_H", "OP_E"),
  "FS" = c("FS_A", "FS_Z"),
  "PF" = c("PF_C", "PF_S"),
  "FT" = c("FT_D", "FT_F"),
  "FD" = c("FD_D", "FD_I")
)

# Validate Input
if (!TRAIT_TO_ANALYZE %in% names(TRAIT_DEFINITIONS)) {
  stop("Error: Invalid trait abbreviation. Please check the list.")
}

# Get current trait levels
current_levels <- TRAIT_DEFINITIONS[[TRAIT_TO_ANALYZE]]

# 4. Path Setup
#-----------------------------------------------------------------------------
tree_file_path <- "GBOTB_Flora_GlobalID.tre"
data_file_path <- "Trait_Data_No_Imputation_QC_Ref.csv"

# Create Output Directory
results_dir <- file.path("Flora_Evolution_Models", TRAIT_TO_ANALYZE)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cat("============================================================\n")
cat("Starting Evolutionary Model Analysis\n")
cat("  - Trait:", TRAIT_TO_ANALYZE, "\n")
cat("  - Order:", paste(current_levels, collapse = " -> "), "(State 1 -> State 2)\n")
cat("  - Output Dir:", results_dir, "\n")
cat("============================================================\n\n")

# 5. Data Loading and Cleaning
#-----------------------------------------------------------------------------
cat("1. Reading Phylogenetic Tree...\n")
if(!file.exists(tree_file_path)) stop("Tree file not found.")
phy <- read.tree(tree_file_path)

cat("2. Reading Trait Data...\n")
if(!file.exists(data_file_path)) stop("Data file not found.")
raw_data <- read.csv(data_file_path, stringsAsFactors = FALSE, check.names = FALSE)

# Select Columns
if (!TRAIT_TO_ANALYZE %in% colnames(raw_data)) {
  stop(paste("Error: Column", TRAIT_TO_ANALYZE, "not found in data."))
}

trait_data <- raw_data %>%
  select(GlobalID, all_of(TRAIT_TO_ANALYZE)) %>%
  rename(Species = GlobalID, Trait_Value = !!TRAIT_TO_ANALYZE)

# Clean Data: Remove NA/Empty, Filter levels, Set Factors
cat("3. Cleaning Data...\n")
trait_data_clean <- trait_data %>%
  filter(!is.na(Trait_Value) & Trait_Value != "") %>%
  filter(Trait_Value %in% current_levels) %>%
  mutate(Trait_Value = factor(Trait_Value, levels = current_levels))

# Intersect Tree and Data
combined_species <- intersect(phy$tip.label, trait_data_clean$Species)
n_species <- length(combined_species)

if (n_species == 0) {
  stop("Error: No overlapping species between tree and data.")
}

# Prune Tree and Data
phy_pruned <- keep.tip(phy, combined_species)
data_pruned <- trait_data_clean %>% filter(Species %in% combined_species)

# Format for corHMM: [Species, Trait]
data_corhmm <- data_pruned %>% select(Species, Trait_Value)

cat(sprintf("   - Number of Species analyzed: %d\n\n", n_species))

# 6. Model Fitting (ER & ARD)
#-----------------------------------------------------------------------------
model_types <- c("ER", "ARD")
model_fits <- list()

for (mod in model_types) {
  
  rds_path <- file.path(results_dir, paste0(mod, "_fit.rds"))
  
  if (file.exists(rds_path)) {
    cat(sprintf("   - Model exists: %s (Skipping calculation, loading file)\n", mod))
    model_fits[[mod]] <- readRDS(rds_path)
  } else {
    cat(sprintf("   - Fitting model: %s ...\n", mod))
    
    # corHMM input: Trait must be integer (1, 2)
    input_data <- data_corhmm
    input_data$Trait_Value <- as.integer(input_data$Trait_Value)
    
    # Define Rate Matrix Structure
    if (mod == "ER") {
      # ER: All off-diagonal rates are the same parameter (1)
      rate_mat <- matrix(1, nrow=2, ncol=2)
      diag(rate_mat) <- NA
    } else {
      # ARD: 1->2 is param 1, 2->1 is param 2
      rate_mat <- matrix(c(NA, 2, 1, NA), nrow=2, byrow=TRUE)
    }
    
    # Run corHMM
    start_t <- Sys.time()
    fit <- corHMM(
      phy = phy_pruned,
      data = input_data,
      rate.cat = 1,          # No hidden states
      rate.mat = rate_mat,
      model = mod,
      node.states = "marginal",
      root.p = "maddfitz",   
      n.cores = 1,
      upper.bound = 1        # Limit rate upper bound
    )
    end_t <- Sys.time()
    
    # Save Results
    saveRDS(fit, rds_path)
    
    # Save text summary
    sink(file.path(results_dir, paste0(mod, "_summary.txt")))
    print(fit)
    sink()
    
    cat(sprintf("     Done. Duration: %.2f mins\n", as.numeric(difftime(end_t, start_t, units = "mins"))))
    model_fits[[mod]] <- fit
  }
}

# 7. Result Aggregation and Comparison
#-----------------------------------------------------------------------------
cat("\nSummarizing Results...\n")

if (length(model_fits) > 0) {
  
  res_table <- data.frame(
    Trait = TRAIT_TO_ANALYZE,
    Model = names(model_fits),
    LogLik = sapply(model_fits, function(x) x$loglik),
    AICc = sapply(model_fits, function(x) x$AICc),
    # Extract Rates: 1->2 (Forward) and 2->1 (Backward)
    Rate_1_to_2 = sapply(model_fits, function(x) x$solution[1,2]),
    Rate_2_to_1 = sapply(model_fits, function(x) x$solution[2,1]),
    stringsAsFactors = FALSE
  ) %>%
    arrange(AICc) %>%
    mutate(Delta_AICc = AICc - min(AICc)) %>%
    mutate(Weight = exp(-0.5 * Delta_AICc) / sum(exp(-0.5 * Delta_AICc)))
  
  # Format Numbers
  res_table_print <- res_table %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
  
  # Save Comparison Table
  out_csv <- file.path(results_dir, paste0("Model_Comparison_", TRAIT_TO_ANALYZE, ".csv"))
  write.csv(res_table_print, out_csv, row.names = FALSE)
  
  cat("------------------------------------------------------------\n")
  cat("Model Comparison (Sorted by AICc):\n")
  print(res_table_print)
  cat("------------------------------------------------------------\n")
  cat(sprintf("State 1 = %s (Ancestral)\n", current_levels[1]))
  cat(sprintf("State 2 = %s (Derived)\n", current_levels[2]))
  cat(sprintf("Rate_1_to_2: Rate from %s to %s\n", current_levels[1], current_levels[2]))
  cat("------------------------------------------------------------\n")
  cat("Analysis Completed.\n")
  
} else {
  cat("Error: No models in list.\n")
}

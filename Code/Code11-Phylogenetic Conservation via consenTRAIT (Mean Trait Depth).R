#-----------------------------------------------------------------------------
# Code 11: Phylogenetic Conservation via consenTRAIT (Mean Trait Depth)
#
# Description:
# Calculates the "Mean Trait Depth" using `castor::consentrait_depth`.
# This metric estimates the average depth (age) of clades sharing a specific trait.
# It is much faster than Fritz's D and suitable for very large trees (60k+ tips).
#
# Inputs:
# - GBOTB_Flora_GlobalID.tre (Tree)
# - Trait_Data_No_Imputation_QC_Ref.csv (Raw Trait Data)
#
# Outputs:
# - Phylogenetic_consenTRAIT.csv
#-----------------------------------------------------------------------------

# 0. Global Settings & Libraries
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("castor", "ape", "dplyr", "data.table")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(castor)
library(ape)
library(dplyr)
library(data.table)

# 1. Parameters
#-----------------------------------------------------------------------------
# Input Files
TREE_PATH <- "GBOTB_Flora_GlobalID.tre"
DATA_PATH <- "Trait_Data_No_Imputation_QC_Ref.csv"

# Output File
OUTPUT_PATH <- "Phylogenetic_consenTRAIT.csv"

# Permutations for P-value calculation
N_PERMUTATIONS <- 1000

# Trait Definitions: c(Background_State, Target_State)
# consenTRAIT calculates the depth of the *second* state (Target_State)
TRAIT_DEFINITIONS <- list(
  "GF" = c("GF_W", "GF_H"), # Target: Herbaceous (H)
  "LC" = c("LC_P", "LC_A"), # Target: Annual (A)
  "PH" = c("PH_E", "PH_D"), # Target: Deciduous (D)
  "RS" = c("RS_T", "RS_F"), # Target: Fibrous root (F)
  "PT" = c("PT_U", "PT_M"), # Target: Multiple phyllotaxy (M)
  "LS" = c("LS_S", "LS_C"), # Target: Compound leaf (C)
  "LM" = c("LM_E", "LM_L"), # Target: Lobed (L)
  "LV" = c("LV_R", "LV_P"), # Target: Parallel veins (P)
  "SS" = c("SS_B", "SS_U"), # Target: Unisexual (U)
  "IF" = c("IF_D", "IF_I"), # Target: Indeterminate (I)
  "OP" = c("OP_H", "OP_E"), # Target: Epigynous (E)
  "FS" = c("FS_A", "FS_Z"), # Target: Zygomorphic (Z)
  "PF" = c("PF_C", "PF_S"), # Target: Sympetalous (S)
  "FT" = c("FT_D", "FT_F"), # Target: Fleshy fruit (F)
  "FD" = c("FD_D", "FD_I")  # Target: Indehiscent (I)
)

cat("========= Starting consenTRAIT Analysis =========\n")

# 2. Data Loading
#-----------------------------------------------------------------------------
cat("1. Loading Phylogenetic Tree...\n")
if(!file.exists(TREE_PATH)) stop("Tree file not found.")
full_tree <- ape::read.tree(TREE_PATH)
full_tree$node.label <- NULL # Clean labels
cat(sprintf("   -> Loaded tree with %d tips.\n", length(full_tree$tip.label)))

cat("2. Loading Trait Data...\n")
if(!file.exists(DATA_PATH)) stop("Data file not found.")
trait_data_raw <- fread(DATA_PATH)

# 3. Main Analysis Loop
#-----------------------------------------------------------------------------
results_list <- list()

cat("\n--- Analyzing 15 Traits ---\n")

for (trait_name in names(TRAIT_DEFINITIONS)) {
  
  cat(sprintf("  -> Processing: %s ... ", trait_name))
  
  # A. Preprocessing
  state_0_label <- TRAIT_DEFINITIONS[[trait_name]][1]
  state_1_label <- TRAIT_DEFINITIONS[[trait_name]][2]
  
  # Filter Data
  # Using data.table syntax for speed
  trait_subset <- trait_data_raw[
    !is.na(get(trait_name)) & get(trait_name) %in% c(state_0_label, state_1_label),
    .(GlobalID, Trait_Val = get(trait_name))
  ]
  
  # Convert to Numeric Binary
  trait_subset[, tip_state_numeric := ifelse(Trait_Val == state_1_label, 1, 0)]
  
  # B. Tree Matching
  common_tips <- intersect(full_tree$tip.label, trait_subset$GlobalID)
  
  if(length(common_tips) < 50) {
    cat("Skipped (Too few species < 50)\n")
    next
  }
  
  # Prune Tree
  analysis_tree <- ape::keep.tip(full_tree, common_tips)
  
  # Align Data to Tree Tip Order
  # Re-order trait_subset based on tree tip labels
  analysis_data <- trait_subset[match(analysis_tree$tip.label, GlobalID)]
  
  # C. Run consenTRAIT
  tip_states_vector <- analysis_data$tip_state_numeric
  
  # Parameters:
  # min_fraction = 0.9: A clade is considered "positive" if 90% of tips are state 1.
  consentrait_results <- tryCatch({
    castor::consentrait_depth(
      tree = analysis_tree,
      tip_states = tip_states_vector,
      min_fraction = 0.9, 
      count_singletons = TRUE, 
      Npermutations = N_PERMUTATIONS
    )
  }, error = function(e) {
    return(NULL)
  })
  
  # D. Store Results
  if (!is.null(consentrait_results)) {
    
    # Interpretation: 
    # High Mean_Depth = Traits fixed deep in the phylogeny (Old/Conserved)
    # Low Mean_Depth = Traits fixed recently (Young/Labile)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Trait = trait_name,
      Target_State = state_1_label,
      NTips = length(analysis_tree$tip.label),
      Mean_Depth = round(consentrait_results$mean_depth, 4),
      P_Value = consentrait_results$P,
      Mean_Random_Depth = round(consentrait_results$mean_random_depth, 4),
      N_Positive_Clades = consentrait_results$Npositives
    )
    cat(sprintf("Done. Depth=%.4f (P=%.3f)\n", 
                consentrait_results$mean_depth, consentrait_results$P))
  } else {
    cat("Failed (Calculation Error)\n")
  }
}

# 4. Export Results
#-----------------------------------------------------------------------------
cat("\n--- Saving Results ---\n")

if (length(results_list) > 0) {
  final_summary_df <- bind_rows(results_list)
  
  # Sort by Depth descending
  final_summary_df <- final_summary_df %>% arrange(desc(Mean_Depth))
  
  fwrite(final_summary_df, OUTPUT_PATH)
  
  cat("Success! Results saved to:\n", OUTPUT_PATH, "\n")
  print(head(final_summary_df))
  
} else {
  cat("\nWarning: No valid results generated.\n")
}

#-----------------------------------------------------------------------------
# Code 10: Phylogenetic Signal of Binary Traits (Fritz's D via Subsampling)
#
# Description:
# Calculates Fritz's D for binary traits using a bootstrapping approach.
# Due to the large size of the phylogenetic tree (~60k+ tips), we use 
# repeated subsampling to estimate the signal.
#
# Strategy:
# 1. For each trait, filter valid species.
# 2. Randomly subsample 1,000 species.
# 3. Calculate D statistic (100 permutations per run).
# 4. Repeat 100 times and calculate Mean & SD.
#
# Inputs:
# - GBOTB_Flora_GlobalID.tre (Phylogenetic Tree)
# - Trait_Data_No_Imputation_QC_Ref.csv (Raw Trait Data)
#
# Outputs:
# - Phylogenetic_Signal_Summary.csv
#-----------------------------------------------------------------------------

# 0. Global Settings & Libraries
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("caper", "dplyr", "ape", "data.table")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(caper)
library(dplyr)
library(ape)
library(data.table)

# 1. Parameters
#-----------------------------------------------------------------------------
# [Bootstrap Settings]
SUBSAMPLE_SIZE <- 1000    # Species count per iteration
N_ITERATIONS   <- 100     # Number of bootstrap iterations
PERMUTATIONS   <- 100     # Permutations inside phylo.d for P-value

# [Trait Definitions]
# Format: Trait_Name = c("State_0", "State_1")
TRAIT_DEFINITIONS <- list(
  "GF" = c("GF_W", "GF_H"), "LC" = c("LC_P", "LC_A"), "PH" = c("PH_E", "PH_D"),
  "RS" = c("RS_T", "RS_F"), "PT" = c("PT_U", "PT_M"), "LS" = c("LS_S", "LS_C"),
  "LM" = c("LM_E", "LM_L"), "LV" = c("LV_R", "LV_P"), "SS" = c("SS_B", "SS_U"),
  "IF" = c("IF_D", "IF_I"), "OP" = c("OP_H", "OP_E"), "FS" = c("FS_A", "FS_Z"),
  "PF" = c("PF_C", "PF_S"), "FT" = c("FT_D", "FT_F"), "FD" = c("FD_D", "FD_I")
)

target_traits <- names(TRAIT_DEFINITIONS)

# 2. Load Data
#-----------------------------------------------------------------------------
tree_file <- "GBOTB_Flora_GlobalID.tre"
data_file <- "Trait_Data_No_Imputation_QC_Ref.csv"

if(!file.exists(tree_file) || !file.exists(data_file)) {
  stop("Error: Input files (Tree or CSV) not found!")
}

cat("1. Reading Phylogenetic Tree...\n")
phy_full <- read.tree(tree_file)
phy_full$node.label <- NULL # Remove node labels to prevent conflict

cat("2. Reading Trait Data...\n")
raw_data <- fread(data_file) 

# Container for final results
final_summary <- data.frame()

# 3. Main Calculation Loop
#-----------------------------------------------------------------------------
cat(sprintf("3. Calculating Phylogenetic Signal for %d traits...\n", length(target_traits)))
cat(sprintf("   [Settings] Sample: %d | Iterations: %d | Permutations: %d\n", 
            SUBSAMPLE_SIZE, N_ITERATIONS, PERMUTATIONS))

for (tr in target_traits) {
  
  cat(sprintf("\n>> Processing: [%s] ... ", tr))
  
  # 3.1 Prepare Data for Current Trait
  current_states <- TRAIT_DEFINITIONS[[tr]] 
  
  # Filter data: Non-NA and matches defined states
  sub_data_full <- raw_data[
    !is.na(get(tr)) & get(tr) != "" & get(tr) %in% current_states, 
    .(GlobalID, Trait_Value = get(tr))
  ]
  
  # Find intersection with tree
  common_sp_full <- intersect(phy_full$tip.label, sub_data_full$GlobalID)
  total_sp_count <- length(common_sp_full)
  
  if(total_sp_count < 20) {
    cat("Skipped (Too few species)\n")
    next
  }
  
  cat(sprintf("(Pool N=%d) -> Bootstrapping:\n", total_sp_count))
  
  # 3.2 Iteration Loop
  iter_results <- data.frame()
  pb <- txtProgressBar(min = 0, max = N_ITERATIONS, style = 3)
  
  for (i in 1:N_ITERATIONS) {
    tryCatch({
      # A. Sample IDs
      if (total_sp_count <= SUBSAMPLE_SIZE) {
        sampled_sp <- common_sp_full
      } else {
        sampled_sp <- sample(common_sp_full, SUBSAMPLE_SIZE)
      }
      
      # B. Prune Tree & Data
      phy_sub <- keep.tip(phy_full, sampled_sp)
      data_sub <- sub_data_full[GlobalID %in% sampled_sp]
      
      # C. Convert to Binary (0/1)
      # State 1 = 0, State 2 = 1
      data_sub[, Binary_Value := ifelse(Trait_Value == current_states[1], 0, 1)]
      
      # D. Check Monomorphism
      if (length(unique(data_sub$Binary_Value)) < 2) next 
      
      # E. Calculate Fritz's D
      # Note: data.frame is required for comparative.data (tibbles/DT can cause issues)
      comp_data <- comparative.data(
        phy = phy_sub,
        data = as.data.frame(data_sub), 
        names.col = "GlobalID",
        vcv = FALSE,
        warn.dropped = FALSE,
        force.root = TRUE
      )
      
      d_stat <- phylo.d(
        data = comp_data,
        binvar = Binary_Value,
        permut = PERMUTATIONS
      )
      
      # Store
      iter_results <- rbind(iter_results, data.frame(
        D_Est = d_stat$DEstimate,
        P_Rand = d_stat$Pval1, # Prob of Evolving Randomly (Brownian Motion)
        P_BM = d_stat$Pval0    # Prob of Brownian Motion
      ))
      
    }, error = function(e) {
      # Fail silently for single iteration errors
    })
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 3.3 Aggregate Results
  if (nrow(iter_results) > 0) {
    
    mean_D <- mean(iter_results$D_Est, na.rm = TRUE)
    sd_D   <- sd(iter_results$D_Est, na.rm = TRUE)
    mean_P_Rand <- mean(iter_results$P_Rand, na.rm = TRUE)
    
    # Interpretation Logic
    interp <- "Unknown"
    if (mean_D <= 0) {
      interp <- "Highly Conserved"
    } else if (mean_D > 0 & mean_D < 0.5) {
      interp <- "Conserved"
    } else if (mean_D >= 0.5 & mean_D < 1) {
      interp <- "Weak Signal"
    } else {
      interp <- "Random/Overdispersed"
    }
    
    final_summary <- rbind(final_summary, data.frame(
      Trait = tr,
      Total_Pool_Size = total_sp_count,
      Sampled_Size = ifelse(total_sp_count < SUBSAMPLE_SIZE, total_sp_count, SUBSAMPLE_SIZE),
      Mean_D_Value = round(mean_D, 4),
      SD_D_Value = round(sd_D, 4),
      Mean_P_Random = round(mean_P_Rand, 4),
      Interpretation = interp,
      State_0_Ref = current_states[1],
      State_1_Alt = current_states[2]
    ))
    
    cat(sprintf("   Done! Mean D=%.3f (+/-%.3f) -> %s\n", mean_D, sd_D, interp))
    
  } else {
    cat("   Failed: No valid iterations.\n")
  }
}

# 4. Save Results
#-----------------------------------------------------------------------------
output_file <- "Phylogenetic_Signal_Summary.csv"

if (nrow(final_summary) > 0) {
  fwrite(final_summary, output_file)
  
  cat("\n============================================================\n")
  cat("Calculation Completed. Results saved to:\n")
  cat(file.path(getwd(), output_file), "\n")
  cat("------------------------------------------------------------\n")
  print(final_summary[, c("Trait", "Mean_D_Value", "SD_D_Value", "Interpretation")])
  cat("============================================================\n")
} else {
  cat("\nNo results generated.\n")
}

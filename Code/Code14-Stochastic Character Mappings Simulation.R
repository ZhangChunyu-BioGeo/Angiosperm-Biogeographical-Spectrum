#-----------------------------------------------------------------------------
# Code 14: Stochastic Character Mappings Simulation
#
# Description:
# 1. Reads 'Evolution_Models_Summary.csv' to get the best model and rates for a specific trait.
# 2. Constructs the Q matrix.
# 3. Performs 1000 Stochastic Character Mappings (SCM).
# 4. Calculates statistics in dynamic time bins (1 Myr).
#
# Inputs:
#   - Evolution_Models_Summary.csv (Aggregated summary of models)
#   - GBOTB_Flora_GlobalID.tre
#   - Trait_Data_Filtered_No_Imputation.csv
#
# Outputs:
#   - Flora_SCM_Simulations/{Trait}/Simulations_{Model}_1Ma/...
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

library(ape)
library(dplyr)
library(phytools)
library(R.utils)

# 2. User Control Panel
#-----------------------------------------------------------------------------

# [Target Trait Abbreviation]
# Options: "GF", "LC", "PH", "RS", "PT", "LS", "LM", "LV", 
#          "SS", "IF", "OP", "FS", "PF", "FT", "FD"
TRAIT_ABBR <- "GF"

# Simulation Parameters
START_SEED <- 1
END_SEED   <- 1000

TIMEOUT_SECONDS  <- 9000
TIME_BIN_SIZE_MA <- 1

# File Names
tree_file_name    <- "GBOTB_Flora_GlobalID.tre"
data_file_name    <- "Trait_Data_Filtered_No_Imputation.csv"
summary_file_name <- "Evolution_Models_Summary.csv"

# 3. Parameter Extraction
#-----------------------------------------------------------------------------
if (!file.exists(summary_file_name)) {
  stop("Error: Summary file 'Evolution_Models_Summary.csv' not found.")
}

# Read Summary
model_summary <- read.csv(summary_file_name, stringsAsFactors = FALSE)

# Filter for Target Trait
target_params <- model_summary %>% filter(Trait == TRAIT_ABBR)

if (nrow(target_params) == 0) {
  stop(paste("Error: Trait", TRAIT_ABBR, "not found in summary file."))
}

# Extract Best Model and Rates
BEST_MODEL       <- target_params$Model[1]       # "ER" or "ARD"
AUTO_RATE_1_TO_2 <- target_params$Rate_1_to_2[1]
AUTO_RATE_2_TO_1 <- target_params$Rate_2_to_1[1]

cat("============================================================\n")
cat("SCM Automation Configuration\n")
cat("------------------------------------------------------------\n")
cat("  - Target Trait:", TRAIT_ABBR, "\n")
cat("  - Best Model:", BEST_MODEL, "\n")
cat("  - Rate (1->2):", AUTO_RATE_1_TO_2, "\n")
cat("  - Rate (2->1):", AUTO_RATE_2_TO_1, "\n")
cat("============================================================\n")

# 4. Trait Definitions (State Order)
#-----------------------------------------------------------------------------
TRAIT_DEFINITIONS <- list(
  "GF" = c("GF_W", "GF_H"), "LC" = c("LC_P", "LC_A"), "PH" = c("PH_E", "PH_D"),
  "RS" = c("RS_T", "RS_F"), "PT" = c("PT_U", "PT_M"), "LS" = c("LS_S", "LS_C"),
  "LM" = c("LM_E", "LM_L"), "LV" = c("LV_R", "LV_P"), "SS" = c("SS_B", "SS_U"),
  "IF" = c("IF_D", "IF_I"), "OP" = c("OP_H", "OP_E"), "FS" = c("FS_A", "FS_Z"),
  "PF" = c("PF_C", "PF_S"), "FT" = c("FT_D", "FT_F"), "FD" = c("FD_D", "FD_I")
)
state_names <- TRAIT_DEFINITIONS[[TRAIT_ABBR]]

# 5. Construct Q Matrix & Setup Paths
#-----------------------------------------------------------------------------

# Output Directories
output_dir <- file.path("Flora_SCM_Simulations", TRAIT_ABBR)
sim_results_dir <- file.path(output_dir, paste0("Simulations_", BEST_MODEL, "_", TIME_BIN_SIZE_MA, "Ma"))
failed_dir <- file.path(output_dir, "Failed_Simulations")

if (!dir.exists(sim_results_dir)) dir.create(sim_results_dir, recursive = TRUE)
if (!dir.exists(failed_dir)) dir.create(failed_dir, recursive = TRUE)

# Construct Q Matrix
Q_matrix <- matrix(0, nrow = 2, ncol = 2)
colnames(Q_matrix) <- rownames(Q_matrix) <- state_names

# Fill Rates
Q_matrix[state_names[1], state_names[2]] <- AUTO_RATE_1_TO_2 # State 1 -> 2
Q_matrix[state_names[2], state_names[1]] <- AUTO_RATE_2_TO_1 # State 2 -> 1
diag(Q_matrix) <- -rowSums(Q_matrix)

cat("Constructed Q Matrix:\n")
print(Q_matrix)
cat("\n")

# 6. Data Loading and Cleaning
#-----------------------------------------------------------------------------
cat("Preparing Data...\n")
phy <- read.tree(tree_file_name)
raw_data <- read.csv(data_file_name, stringsAsFactors = FALSE, check.names = FALSE)

# Filter and Clean
trait_data <- raw_data %>%
  dplyr::select(GlobalID, dplyr::all_of(TRAIT_ABBR)) %>%
  dplyr::rename(Species = GlobalID, Trait_Value = !!TRAIT_ABBR) %>%
  dplyr::filter(!is.na(Trait_Value) & Trait_Value != "") %>%
  dplyr::filter(Trait_Value %in% state_names)

valid_species <- intersect(phy$tip.label, trait_data$Species)
if (length(valid_species) == 0) stop("Error: No overlapping species found!")

phy_pruned <- keep.tip(phy, valid_species)
data_pruned <- trait_data %>% dplyr::filter(Species %in% valid_species)

# Convert to Named Vector for make.simmap
x_vector <- setNames(data_pruned$Trait_Value, data_pruned$Species)
x_vector <- x_vector[phy_pruned$tip.label]

# 7. Dynamic Binning Setup
#-----------------------------------------------------------------------------
tree_height <- max(node.depth.edgelength(phy_pruned))
time_break_points <- seq(0, ceiling(tree_height), by = TIME_BIN_SIZE_MA)
num_bins <- length(time_break_points) - 1
node_ages <- tree_height - node.depth.edgelength(phy_pruned)

# 8. Helper Function
#-----------------------------------------------------------------------------
create_empty_stats <- function() {
  df <- data.frame(
    bin_id = 1:num_bins,
    bin_start_age = time_break_points[-(num_bins + 1)],
    bin_end_age = time_break_points[-1],
    total_branch_length = 0
  )
  for (s in state_names) df[[paste0("len_", s)]] <- 0
  df[[paste0("trans_", state_names[1], "_to_", state_names[2])]] <- 0
  df[[paste0("trans_", state_names[2], "_to_", state_names[1])]] <- 0
  return(df)
}

# 9. SCM Simulation Main Loop
#-----------------------------------------------------------------------------
cat("Starting SCM Simulations (", START_SEED, " - ", END_SEED, ")...\n")

for (seed in START_SEED:END_SEED) {
  
  file_name <- sprintf("scm_sim_%s_%s_seed_%04d.csv", TRAIT_ABBR, BEST_MODEL, seed)
  out_path <- file.path(sim_results_dir, file_name)
  
  if (file.exists(out_path)) {
    if (seed %% 100 == 0) cat(sprintf("Seed %d exists, skipping...\n", seed))
    next
  }
  if (seed %% 10 == 0) cat(sprintf("Seed %d ...\n", seed))
  
  # --- Run Simulation ---
  set.seed(seed)
  sim_res <- NULL
  
  tryCatch({
    sim_res <- withTimeout({
      phytools::make.simmap(
        tree = phy_pruned, 
        x = x_vector, 
        Q = Q_matrix, 
        pi = "estimated", 
        nsim = 1, 
        message = FALSE
      )
    }, timeout = TIMEOUT_SECONDS, onTimeout = "error")
  }, error = function(e) {
    err_msg <- gsub("\n", " ", e$message)
    cat(sprintf("Seed %d Failed: %s\n", seed, err_msg))
    write.table(data.frame(Seed=seed, Error=err_msg), file.path(failed_dir, "errors.csv"), 
                append=TRUE, sep=",", row.names=FALSE, col.names=!file.exists(file.path(failed_dir, "errors.csv")))
  })
  
  if (is.null(sim_res)) next
  
  # --- Binning Statistics ---
  current_stats <- create_empty_stats()
  edges <- sim_res$edge
  maps <- sim_res$maps
  
  for (i in 1:nrow(edges)) {
    parent_node <- edges[i, 1]
    start_age <- node_ages[parent_node] 
    edge_map <- maps[[i]]
    current_seg_start <- start_age
    
    for (j in 1:length(edge_map)) {
      state <- names(edge_map)[j]
      duration <- edge_map[j]
      current_seg_end <- current_seg_start - duration
      
      overlap_indices <- which(current_stats$bin_start_age < current_seg_start & 
                                 current_stats$bin_end_age > current_seg_end)
      
      for (bin_idx in overlap_indices) {
        ov_start <- max(current_stats$bin_start_age[bin_idx], current_seg_end)
        ov_end   <- min(current_stats$bin_end_age[bin_idx], current_seg_start)
        ov_len   <- ov_end - ov_start
        
        current_stats$total_branch_length[bin_idx] <- current_stats$total_branch_length[bin_idx] + ov_len
        col_name <- paste0("len_", state)
        current_stats[bin_idx, col_name] <- current_stats[bin_idx, col_name] + ov_len
      }
      
      if (j > 1) {
        prev_state <- names(edge_map)[j-1]
        trans_time <- current_seg_start
        trans_bin_idx <- which(current_stats$bin_start_age <= trans_time & 
                                 current_stats$bin_end_age > trans_time)
        if (length(trans_bin_idx) > 0) {
          idx <- trans_bin_idx[1]
          trans_col <- paste0("trans_", prev_state, "_to_", state)
          current_stats[idx, trans_col] <- current_stats[idx, trans_col] + 1
        }
      }
      current_seg_start <- current_seg_end
    }
  }
  
  # --- Calculate Proportions and Rates ---
  s1 <- state_names[1]
  s2 <- state_names[2]
  
  current_stats[[paste0("prop_", s1)]] <- current_stats[[paste0("len_", s1)]] / current_stats$total_branch_length
  current_stats[[paste0("prop_", s2)]] <- current_stats[[paste0("len_", s2)]] / current_stats$total_branch_length
  
  # Rate Calculation (Transitions per Unit Branch Length)
  len_s1 <- current_stats[[paste0("len_", s1)]]
  trans_12 <- current_stats[[paste0("trans_", s1, "_to_", s2)]]
  current_stats[[paste0("rate_", s1, "_to_", s2)]] <- ifelse(len_s1 > 0, trans_12 / len_s1, 0)
  
  len_s2 <- current_stats[[paste0("len_", s2)]]
  trans_21 <- current_stats[[paste0("trans_", s2, "_to_", s1)]]
  current_stats[[paste0("rate_", s2, "_to_", s1)]] <- ifelse(len_s2 > 0, trans_21 / len_s2, 0)
  
  write.csv(current_stats, out_path, row.names = FALSE)
}

cat("\n=======================================================\n")
cat("SCM Simulation Task Completed.\n")
cat("Output Directory:", sim_results_dir, "\n")
cat("=======================================================\n")

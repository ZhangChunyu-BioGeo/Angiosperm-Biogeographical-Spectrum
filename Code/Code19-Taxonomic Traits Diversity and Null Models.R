#-----------------------------------------------------------------------------
# Code 19: Taxonomic Traits Diversity and Null Models
#
# Note on Data Source:
#   Species distribution data is derived from the Global Biodiversity Information 
#   Facility (GBIF) (https://www.gbif.org/). 
#   This repository provides a demonstration dataset: 
#   "[DEMO]Global_Ecoregion_Occ_Plants_GBIF.csv", which is a subsample of 
#   100,000 records.
#   The download record for this study is available at GBIF: https://doi.org/10.15468/dl.dj8zxb.
#
# IMPORTANT: 
#   To fully reproduce this analysis, please follow the instructions to download 
#   species distribution data from GBIF and aggregate it to ecoregions. 
#   All other data required for the main analysis after species aggregation 
#   is provided in this repository.
#
#-----------------------------------------------------------------------------

# 1. Environment and Parameter Settings
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

# --- Input Filename Setting ---
# Change this to the full filename when running the complete analysis
OCC_FILENAME <- "[DEMO]Global_Ecoregion_Occ_Plants_GBIF.csv" 

# --- Core Control Parameters ---
N_PERM      <- 1000  # Number of Null Model Permutations

# --- Task Batching Settings ---
# Divide total tasks into N_PARTS
N_PARTS     <- 10   
# Set Start and End parts for this specific run (e.g., 1-10 for all)
START_PART  <- 1
END_PART    <- 10

# --- Automatic DEMO Detection & Output Setup ---
is_demo_mode <- grepl("^\\[DEMO\\]", OCC_FILENAME)

# Define Base Output Directory
output_dir_name <- "Taxonomic_Traits_Diversity_Null_Model"

if (is_demo_mode) {
  # Trigger Warning as requested
  warning("----------------------------------------------------------------\n",
          " [WARNING] RUNNING IN DEMO MODE \n",
          " Detected input: ", OCC_FILENAME, "\n",
          " Output directory and files will be prefixed with [DEMO].\n",
          " Please use the full GBIF dataset for final analysis.\n",
          "----------------------------------------------------------------")
  
  # Update Output Directory Name
  output_dir_name <- paste0("[DEMO]", output_dir_name)
}

# Create Output Directory
if (!dir.exists(output_dir_name)) dir.create(output_dir_name)

# Load Packages
pkgs <- c("data.table", "dplyr", "fastmatch")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(dplyr)
library(fastmatch)

# 2. Data Loading and Pre-processing
#-----------------------------------------------------------------------------
cat(">>> 1. Loading and Pre-processing Data...\n")

# A. Trait Data (No Imputation, QC Reference)
# ----------------------------------------------------------------
trait_path <- "Trait_Data_No_Imputation_QC_Ref.csv"
if(!file.exists(trait_path)) stop("Trait data not found.")

trait_df <- fread(trait_path)
trait_cols <- c("GF", "LC", "PH", "RS", "PT", "LS", "LM", "LV", 
                "SS", "IF", "OP", "FS", "PF", "FT", "FD")

# Convert to 0/1 Numeric Matrix (Core Optimization)
cat("   -> Converting traits to binary numeric matrix...\n")
trait_mat <- trait_df[, ..trait_cols] %>%
  mutate(across(everything(), ~as.numeric(as.factor(.)) - 1)) %>%
  as.matrix()
rownames(trait_mat) <- trait_df$GlobalID

# Ensure no NAs
if (any(is.na(trait_mat))) stop("Error: NA detected in trait matrix!")

# Calculate Global Centroid (for TaxPos)
global_centroid <- colMeans(trait_mat)

# B. Distribution Data
# ----------------------------------------------------------------
if(!file.exists(OCC_FILENAME)) stop(paste("Occurrence file not found:", OCC_FILENAME))

occ_df <- fread(OCC_FILENAME) 
occ_df <- occ_df[, .(GlobalID, ECO_ID)]

# C. Ecoregion Info (Realm)
# ----------------------------------------------------------------
coord_df <- fread("Global_Ecoregion_Coordinates.csv")
eco_info <- coord_df[, .(ECO_ID, Realm)]

# Merge Info
occ_merged <- merge(occ_df, eco_info, by = "ECO_ID", all.x = TRUE)
occ_merged <- occ_merged[!is.na(Realm)] # Remove records with no Realm info

# D. Build Realm-Specific Species Pools
# ----------------------------------------------------------------
cat("   -> Building Realm-specific species pools...\n")
realm_pools <- split(occ_merged$GlobalID, occ_merged$Realm)
realm_pools <- lapply(realm_pools, unique)

# E. Prepare Ecoregion List
# ----------------------------------------------------------------
eco_stats <- occ_merged[, .N, by = .(ECO_ID, Realm)]
full_eco_list <- eco_stats[N >= 5]$ECO_ID # Filter ecoregions with <5 species

# 3. Task Allocation Logic (Batching)
#-----------------------------------------------------------------------------
total_tasks <- length(full_eco_list)
cat(sprintf("   -> Total valid ecoregions detected: %d\n", total_tasks))

# Create Partition Indices
part_indices <- cut(seq_along(full_eco_list), breaks = N_PARTS, labels = FALSE)

# Filter for current run
target_indices <- which(part_indices >= START_PART & part_indices <= END_PART)
run_list <- full_eco_list[target_indices]

cat("========================================================\n")
cat(sprintf("   Task Configuration: Running Parts %d to %d (of %d)\n", START_PART, END_PART, N_PARTS))
cat(sprintf("   Ecoregions to process in this run: %d\n", length(run_list)))
cat(sprintf("   Output Directory: %s\n", output_dir_name))
cat("========================================================\n")

if(length(run_list) == 0) stop("No ecoregions selected in this batch range!")

# 4. Define Core Calculation Function
#-----------------------------------------------------------------------------
calc_tax_metrics <- function(spp_ids, trait_matrix, g_centroid) {
  
  valid_ids <- intersect(spp_ids, rownames(trait_matrix))
  if (length(valid_ids) < 2) return(rep(NA, 4))
  
  sub_mat <- trait_matrix[valid_ids, , drop = FALSE]
  n_spp <- nrow(sub_mat)
  n_traits <- ncol(sub_mat)
  
  # --- 1. TaxPos & TaxDis ---
  local_centroid <- colMeans(sub_mat)
  # TaxPos: Distance from local centroid to global centroid
  tax_pos <- sqrt(sum((local_centroid - g_centroid)^2))
  
  # TaxDis: Mean distance of species to local centroid
  centered_mat <- sweep(sub_mat, 2, local_centroid, "-")
  dists_to_cen <- sqrt(rowSums(centered_mat^2))
  tax_dis <- mean(dists_to_cen)
  
  # --- 2. TaxMPD & TaxMNTD (Hamming/Manhattan for Binary) ---
  d_dist <- dist(sub_mat, method = "manhattan")
  
  # MPD: Mean Pairwise Distance (Normalized by number of traits)
  tax_mpd <- mean(d_dist) / n_traits 
  
  # MNTD: Mean Nearest Taxon Distance
  d_mat <- as.matrix(d_dist) / n_traits
  diag(d_mat) <- Inf
  min_dists <- apply(d_mat, 1, min)
  tax_mntd <- mean(min_dists)
  
  return(c(TaxMPD = tax_mpd, TaxMNTD = tax_mntd, TaxDis = tax_dis, TaxPos = tax_pos))
}

# 5. Main Loop (Serial with Skip Check)
#-----------------------------------------------------------------------------
cat(">>> 2. Starting Serial Processing...\n")

pb <- txtProgressBar(min = 0, max = length(run_list), style = 3)

for (i in 1:length(run_list)) {
  
  curr_eco <- run_list[i]
  setTxtProgressBar(pb, i)
  
  # --- Check Existing (Skip Logic) ---
  sub_dir <- file.path(output_dir_name, curr_eco)
  summary_file <- file.path(sub_dir, "summary_metrics.csv")
  
  if (file.exists(summary_file)) {
    # Skip if file exists
    next
  }
  
  # --- Calculation ---
  if (!dir.exists(sub_dir)) dir.create(sub_dir)
  
  curr_data <- occ_merged[ECO_ID == curr_eco]
  curr_spp <- curr_data$GlobalID
  curr_realm <- curr_data$Realm[1]
  richness <- length(curr_spp)
  
  curr_pool <- realm_pools[[curr_realm]]
  
  # Skip if pool is smaller than richness (rare edge case)
  if (length(curr_pool) < richness) next
  
  # A. Observed Metrics
  obs <- calc_tax_metrics(curr_spp, trait_mat, global_centroid)
  if (is.na(obs[1])) next 
  
  # B. Null Model Simulations
  null_mat <- matrix(NA, nrow = N_PERM, ncol = 4)
  colnames(null_mat) <- c("TaxMPD", "TaxMNTD", "TaxDis", "TaxPos")
  
  for (p in 1:N_PERM) {
    sim_spp <- sample(curr_pool, richness)
    null_mat[p, ] <- calc_tax_metrics(sim_spp, trait_mat, global_centroid)
  }
  
  # C. Save Raw Simulations
  fwrite(as.data.table(null_mat), file.path(sub_dir, "null_simulations.csv"))
  
  # D. Stats & Summary (SES Calculation)
  null_mean <- colMeans(null_mat, na.rm = TRUE)
  null_sd   <- apply(null_mat, 2, sd, na.rm = TRUE)
  
  null_sd[null_sd == 0] <- NA
  ses_val <- (obs - null_mean) / null_sd
  
  res_dt <- data.table(
    ECO_ID = curr_eco,
    Realm = curr_realm,
    Richness = richness,
    
    TaxMPD_Obs = obs["TaxMPD"], TaxMPD_SES = ses_val["TaxMPD"],
    TaxMNTD_Obs = obs["TaxMNTD"], TaxMNTD_SES = ses_val["TaxMNTD"],
    TaxDis_Obs = obs["TaxDis"], TaxDis_SES = ses_val["TaxDis"],
    TaxPos_Obs = obs["TaxPos"], TaxPos_SES = ses_val["TaxPos"],
    
    TaxMPD_Null_Mean = null_mean["TaxMPD"], TaxMPD_Null_SD = null_sd["TaxMPD"],
    TaxMNTD_Null_Mean = null_mean["TaxMNTD"], TaxMNTD_Null_SD = null_sd["TaxMNTD"],
    TaxDis_Null_Mean = null_mean["TaxDis"], TaxDis_Null_SD = null_sd["TaxDis"],
    TaxPos_Null_Mean = null_mean["TaxPos"], TaxPos_Null_SD = null_sd["TaxPos"]
  )
  
  fwrite(res_dt, summary_file)
  
  # Garbage Collection (Prevent Memory Leaks)
  if(i %% 50 == 0) gc()
}

close(pb)

# 6. Aggregation (Combine Results)
#-----------------------------------------------------------------------------
cat("\n>>> 3. Aggregating ALL available results from disk...\n")

# Get all subdirectories in the output root
all_subdirs <- list.dirs(output_dir_name, full.names = TRUE, recursive = FALSE)

# Find valid summary files
valid_files <- file.path(all_subdirs, "summary_metrics.csv")
valid_files <- valid_files[file.exists(valid_files)]

if (length(valid_files) > 0) {
  # Read and bind
  all_res_list <- lapply(valid_files, fread)
  final_df <- rbindlist(all_res_list, fill = TRUE)
  
  # Determine Output Filename
  final_filename <- "Diversity_Taxonomic_Traits_SES_Summary.csv"
  if (is_demo_mode) {
    final_filename <- paste0("[DEMO]", final_filename)
  }
  
  final_output_path <- file.path(output_dir_name, final_filename)
  fwrite(final_df, final_output_path)
  
  cat(sprintf(">>> Summary Updated. Contains %d ecoregions.\n", nrow(final_df)))
  cat(">>> File saved to:", final_output_path, "\n")
  
} else {
  cat(">>> No completed results found yet.\n")
}

cat(">>> Batch Job Finished.\n")

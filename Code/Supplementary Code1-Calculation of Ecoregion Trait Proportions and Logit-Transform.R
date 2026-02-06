#-----------------------------------------------------------------------------
# Supplementary Code 1: Calculation of Ecoregion Trait Proportions and Logit-Transform
#
# [IMPORTANT NOTE]:
# Species distribution data are sourced from the Global Biodiversity Information Facility (GBIF). 
# This repository only provides demonstration distribution data.
# To fully replicate this analysis, please follow the instructions to download 
# species distribution data from GBIF and aggregate it to ecoregions.
# The processed data used for the main analysis of this study are already provided in the repository.
#
# Description:
# This script calculates trait proportions based on species traits and their distribution 
# within ecoregions:
# 1. Cleans distribution data, removing invalid or irrelevant ecoregions.
# 2. Calculates trait proportions and species counts, saving raw results:
#     Trait_Proportion_Ecoregion_No_Imputation.csv
# 3. Filters ecoregions where species counts for all traits are greater than 30, saving results:
#     Trait_Proportion_Ecoregion_No_Imputation_30Greater.csv
# 4. Performs Logit-Transform on filtered proportion data, saving transformed results:
#     Trait_Logit_Proportion_Ecoregion_No_Imputation_30Greater.csv
#
# Data:
# 1. Species distribution data (Global_Ecoregion_Occ_Plants_GBIF.csv) derived from:
#     Global Biodiversity Information Facility (GBIF) (https://www.gbif.org/).
#     The download record for this study is available at: https://doi.org/10.15468/dl.dj8zxb.
#     [DEMO]Global_Ecoregion_Occ_Plants_GBIF.csv is a subsampled demo dataset (100k records).
# 2. Species trait data (refer to README.md or the original article for sources). 
#    Three versions are provided:
#     Trait_Data_No_Imputation_QC_Ref.csv (Main analyses): No imputation.
#     Trait_Data_Imputed_Phylo_Free.csv (Trait divergence): Imputation without phylogeny.
#     Trait_Data_Imputed_Phylo_Informed.csv (Supplementary comparing): Imputation with phylogeny.
# 3. Biome types provided by WWF, merged from Global_Ecoregion_Biome_Type.csv.
#
# Last Update: 2026/02/06
#-----------------------------------------------------------------------------

# 1. Set Working Directory and Load Libraries
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")
library(data.table)

# 2. Parameters and File Settings
#-----------------------------------------------------------------------------
# --- Input Files ---
# Note: Modify this filename if using the full GBIF dataset
occ_matrix_input_file <- "[DEMO]Global_Ecoregion_Occ_Plants_GBIF.csv"

trait_data_input_file <- "Trait_Data_No_Imputation_QC_Ref.csv"
biome_data_input_file <- "Global_Ecoregion_Biome_Type.csv"

# --- Output Files ---
# Automatically determine whether to add [DEMO] prefix based on input filename
ecoregion_trait_output_file    <- "[DEMO]Trait_Proportion_Ecoregion_No_Imputation.csv"
ecoregion_filtered_output_file <- "[DEMO]Trait_Proportion_Ecoregion_No_Imputation_30Greater.csv"
ecoregion_logit_output_file    <- "[DEMO]Trait_Logit_Proportion_Ecoregion_No_Imputation_30Greater.csv"

# --- Target Trait Columns Definition ---
target_trait_cols <- c("GF_H", "LC_A", "PH_D", "RS_T", "PT_U", "LS_C", 
                       "LM_L", "LV_P", "SS_B", "IF_I", "OP_H", "FS_Z", 
                       "PF_C", "FT_D", "FD_D")


# 3. Step 1: Clean Distribution Data
#-----------------------------------------------------------------------------
cat("Step 1/4: Cleaning raw distribution matrix...\n")

if (!file.exists(occ_matrix_input_file)) {
  stop(sprintf("Error: Input file '%s' not found. Check path or filename.", occ_matrix_input_file))
}

eco_to_exclude <- c("-9999", "-9998", "21101")
occ_matrix <- fread(occ_matrix_input_file)
cat(sprintf(" -> Raw matrix contains %s records.\n", format(nrow(occ_matrix), big.mark=",")))

# Filter in memory, do not save intermediate files
occ_matrix_filtered <- occ_matrix[!ECO_ID %in% eco_to_exclude]
cat(sprintf(" -> Remaining records after removing %s ecoregions: %s.\n", length(eco_to_exclude), format(nrow(occ_matrix_filtered), big.mark=",")))


# 4. Step 2: Calculate Trait Proportions and Species Counts per Ecoregion
#-----------------------------------------------------------------------------
cat("Step 2/4: Calculating relative trait proportions and species counts per ecoregion...\n")

# --- 4.1 Load and Prepare Data ---
cat(" -> Loading species trait data...\n")
trait_cols <- c("GlobalID", "GF", "LC", "PH", "RS", "PT", "LS", "LM", "LV", "SS", "IF", "OP", "FS", "PF", "FT", "FD")
trait_data <- fread(trait_data_input_file, select = trait_cols)

cat(" -> Pre-processing trait data...\n")
trait_cols_to_clean <- setdiff(trait_cols, "GlobalID")
for (col in trait_cols_to_clean) {
    prefix <- paste0(col, "_")
    trait_data[, (col) := gsub(prefix, "", get(col), fixed = TRUE)]
}

# --- 4.2 Merge Distribution and Trait Data ---
cat(" -> Merging distribution and trait data...\n")
merged_data <- merge(occ_matrix_filtered, trait_data, by = "GlobalID", all.x = TRUE)
rm(occ_matrix, occ_matrix_filtered, trait_data)
gc()

# --- 4.3 Define Trait Mapping ---
trait_map <- list(
  GF = c("H", "W"), LC = c("A", "P"), PH = c("D", "E"), RS = c("F", "T"),
  PT = c("M", "U"), LS = c("C", "S"), LM = c("E", "L"), LV = c("P", "R"),
  SS = c("B", "U"), IF = c("D", "I"), OP = c("E", "H"), FS = c("A", "Z"),
  PF = c("C", "S"), FT = c("D", "F"), FD = c("D", "I")
)

# --- 4.4 Core Calculation ---
cat(" -> Starting group calculation by ecoregion...\n")
ecoregion_traits <- merged_data[, {
    all_trait_results <- lapply(names(trait_map), function(trait_col_name) {
        valid_traits <- na.omit(.SD[[trait_col_name]])
        n_valid <- length(valid_traits)
        trait_levels <- trait_map[[trait_col_name]]
        
        proportions <- as.list(rep(NA_real_, length(trait_levels)))
        if (n_valid > 0) {
            counts <- table(factor(valid_traits, levels = trait_levels))
            proportions <- as.list(counts / n_valid)
        }
        names(proportions) <- paste(trait_col_name, trait_levels, sep = "_")
        
        species_count <- list(n_valid)
        names(species_count) <- paste0(trait_col_name, "_SpCo")
        
        return(c(proportions, species_count))
    })
    do.call(c, all_trait_results)
}, by = ECO_ID, .SDcols = names(trait_map)]

setorder(ecoregion_traits, ECO_ID)
fwrite(ecoregion_traits, ecoregion_trait_output_file)
cat(sprintf(" -> Raw statistical results saved to: %s\n", ecoregion_trait_output_file))


# 5. Step 3: Filter Data, Merge Biome, and Organize Columns
#-----------------------------------------------------------------------------
cat("Step 3/4: Filtering data (>30 species) and merging Biome information...\n")

# --- 5.1 Filter Species Count > 30 ---
spco_cols <- names(ecoregion_traits)[endsWith(names(ecoregion_traits), "_SpCo")]
all_counts_greater_30 <- Reduce(`&`, lapply(spco_cols, function(col) ecoregion_traits[[col]] > 30))
ecoregion_filtered <- ecoregion_traits[all_counts_greater_30, ]
cat(sprintf(" -> Ecoregions remaining after filtering: %d.\n", nrow(ecoregion_filtered)))

if (nrow(ecoregion_filtered) > 0) {
    # --- 5.2 Calculate SpeciesNumber ---
    ecoregion_filtered[, SpeciesNumber := rowMeans(.SD), .SDcols = spco_cols]
    
    # --- 5.3 Load and Merge Biome ---
    cat(" -> Merging BiomeType...\n")
    if(file.exists(biome_data_input_file)){
        biome_df <- fread(biome_data_input_file)
        # Ensure ECO_ID format consistency
        ecoregion_filtered[, ECO_ID := as.character(ECO_ID)]
        biome_df[, ECO_ID := as.character(ECO_ID)]
        
        ecoregion_filtered <- merge(ecoregion_filtered, biome_df, by = "ECO_ID", all.x = TRUE)
    } else {
        warning(sprintf("Biome file %s not found, skipping merge.", biome_data_input_file))
        ecoregion_filtered[, BiomeType := NA]
    }
    
    # --- 5.4 Extract Specific Columns and Sort ---
    # Define final column order: ECO_ID, SpeciesNumber, BiomeType, then trait columns
    final_cols <- c("ECO_ID", "SpeciesNumber", "BiomeType", target_trait_cols)
    
    # Check for missing columns
    missing_cols <- setdiff(final_cols, names(ecoregion_filtered))
    if(length(missing_cols) > 0) stop(paste("Error: Missing columns", paste(missing_cols, collapse=", ")))
    
    final_data <- ecoregion_filtered[, ..final_cols]
    
    fwrite(final_data, ecoregion_filtered_output_file)
    cat(sprintf(" -> Filtered proportion data saved to: %s\n", ecoregion_filtered_output_file))
    
    
    # 6. Step 4: Logit Transform and Save
    #-----------------------------------------------------------------------------
    cat("Step 4/4: Performing Logit transformation and saving final results...\n")
    
    # --- 6.1 Extract Trait Matrix ---
    trait_mat <- as.matrix(final_data[, ..target_trait_cols])
    
    # --- 6.2 Logit Transformation ---
    cat(" -> Executing Logit transform (delta = 0.001)...\n")
    delta <- 0.001 
    trait_logit <- log((trait_mat + delta) / (1 - trait_mat + delta))
    
    # --- 6.3 Combine Final Data Table ---
    # Keep ECO_ID, SpeciesNumber, BiomeType, replace trait columns with transformed values
    logit_data <- cbind(
        final_data[, .(ECO_ID, SpeciesNumber, BiomeType)],
        as.data.table(trait_logit)
    )
    
    # --- 6.4 Save ---
    fwrite(logit_data, ecoregion_logit_output_file)
    cat(sprintf(" -> Logit transformation results saved to: %s\n", ecoregion_logit_output_file))

} else {
    stop("Error: No ecoregions matched the filter criteria.")
}

cat("\nProcessing complete.\n")

# 7. Check for Demo File and Print Notice
#-----------------------------------------------------------------------------
if (grepl("^\\[DEMO\\]", basename(occ_matrix_input_file))) {
    cat("\n==============================================================\n")
    cat("IMPORTANT NOTE: \n")
    cat("# Species distribution data are sourced from the Global Biodiversity Information Facility (GBIF).\n")
    cat("# This repository only provides demonstration distribution data.\n")
    cat("# To fully replicate this analysis, please follow the instructions to download species\n")
    cat("# distribution data from GBIF and aggregate it to ecoregions.\n")
    cat("# The processed data used for the main analysis of this study are already provided in the repository.\n")
    cat("==============================================================\n")
}

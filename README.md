# README: Code for "Biogeographical Variation of Angiosperm Taxonomic Traits Exhibits a Globally Consistent Structure"

This repository contains the R code and supporting data files for the analyses presented in the manuscript:

> Zhang, C., Zhang, Q., Sun, X., Wang, R., Wang, H. & Zheng, P. (submitted to *Nature Plants*). Biogeographical Variation of Angiosperm Taxonomic Traits Exhibits a Globally Consistent Structure.
> 
> *Hui Wang and Peiming Zheng contributed equally to this study.*

The code implements data processing, principal component analysis (PCA), phylogenetic generalized linear mixed models (PGLMM), phylogenetic signal and conservation analyses (Fritz's D, consenTRAIT), evolutionary model fitting, stochastic character mapping, beta regression, null model simulations for diversity metrics, Generalized Additive Models (GAM) for environmental and trait space mapping, and spatial piecewise structural equation modeling (pSEM).

All analyses were performed in R (version 4.5.1). The code is structured with individual R scripts for each figure in the manuscript, ensuring reproducibility of the global biogeographical spectrum of taxonomic traits across 323,681 angiosperm species and 652 ecoregions. Data files are stored in the `Data` folder.

## License
This project is covered under the MIT License, which allows free use, modification, and distribution. See the LICENSE file for details.
Please note that some data in this repository (e.g., species distribution data from GBIF) may be subject to the licensing terms of their original providers.

## Repository Link
The repository is available on GitHub: https://github.com/ZhangChunyu-BioGeo/Angiosperm-Biogeographical-Spectrum. A static version of this repository is also available on Zenodo: https://doi.org/10.5281/zenodo.20595813.

## System Requirements
- **Operating Systems**: Tested on Windows 10/11, macOS.

- **Software Dependencies**: R (version 4.4.3 or higher). Required R packages (with tested versions): 
ape (5.8), betareg (3.2.3), biscale (1.0.0), caper (1.0.3), castor (1.8.3), corHMM (2.8), cowplot (1.2.0), data.table (1.18.0), dplyr (1.1.4), factoextra (1.0.7), fastmatch (1.1.6), ggplot2 (4.0.0), ggrepel (0.9.6), grid (4.4.1), gridExtra (2.3), hexbin (1.28.5), Matrix (1.7.0), mgcv (1.9.1), nlme (3.1.164), patchwork (1.3.0), pheatmap (1.0.12), phyr (1.1.0), phytools (2.4.4), piecewiseSEM (2.3.1), R.utils (2.12.3), readr (2.1.5), reshape2 (1.4.4), scales (1.4.0), sp (2.1.4), stringr (1.5.1), and tidyr (1.3.1).

- **Hardware**: Standard desktop computer with at least 8 GB RAM and multi-core CPU (e.g., Intel i5 or equivalent). No non-standard hardware required.

- **Tested Versions**: The code has been tested on R 4.4.3 to 4.5.1.

## Installation Guide
1. Install R from https://cran.r-project.org/ and optionally RStudio from https://posit.co/download/rstudio-desktop/.

2. Clone this repository or download and unzip the ZIP file.

3. Open R/RStudio and set the working directory to the repository root.

4. Install dependencies by running:  
   ```R
   install.packages(c(  "data.table", "dplyr", "ggplot2", "ape", "tidyr", "stringr", "phyr",
   "Matrix", "grid", "scales", "factoextra", "gridExtra", "cowplot", "patchwork",
   "phytools", "R.utils", "readr", "mgcv", "sp", "pheatmap", "reshape2",
   "biscale", "ggrepel", "hexbin", "caper", "castor", "corHMM", "betareg",
   "fastmatch", "piecewiseSEM", "nlme"))
   ```

5. Typical installation time: 5-20 minutes on a standard desktop computer with internet access.

## Data Description

1. **Trait Data Sources**: The file `Data/Supplementary_Data/MetaData_of_Synthesis.csv` contains the list of sources for trait data. The `Ref_ID` column in this file corresponds to the `Ref_ID` column in the species trait dataset `Trait_Data_No_Imputation_QC_Ref.csv`.

2. **Species Distribution Data**: The species occurrence records were derived from the **Global Biodiversity Information Facility (GBIF)** (https://www.gbif.org/). This repository only provides a demonstration dataset. The download record for this study is available at GBIF: https://doi.org/10.15468/dl.dj8zxb).

   * The file `[DEMO]Global_Ecoregion_Occ_Plants_GBIF.csv` included in the repository is a subsampled demonstration dataset containing 100,000 records.
   * **Important Note**: To fully reproduce this analysis, please follow the instructions to download the full species distribution data from GBIF and aggregate it to the ecoregion level. All other datasets required for the main analyses (post-species aggregation) are provided in this repository.

3. **Trait Data Versions**: For sources and definitions of species trait data, please refer to the original article. Three versions of the trait dataset are provided:

   * `Trait_Data_No_Imputation_QC_Ref.csv` (**Main analyses**): Non-imputed data (Raw).
   * `Trait_Data_Imputed_Phylo_Free.csv` (**Trait divergence**): Imputed data without phylogenetic information.
   * `Trait_Data_Imputed_Phylo_Informed.csv` (**Supplementary comparing**): Imputed data with phylogenetic information.

4. **Biome Types**: Biome classifications were obtained from the the World Wildlife Fund (WWF). The consolidated data is provided in `Global_Ecoregion_Biome_Type.csv`.

5. **Auxiliary Spatial Data**: Other spatial datasets included are:

   * `Global_Ecoregion_Environment.csv` (Environmental variables)
   * `Global_Ecoregion_Coordinates.csv` (Ecoregion coordinates and area)
   * `Global_Ecoregion_Phylogenetic_Similarity_Matrix_PhyloSor.csv` (Phylogenetic similarity matrix)
   * `Global_Ecoregion_Space_Similarity_Matrix_IDW.csv` (Spatial autocorrelation matrix)
   * These processed datasets are provided in the repository; please refer to the original Article for the specific methodologies used.

6. **Phylogenetic Tree**: The GBOTB tree was derived from the **Open Tree of Life** platform (https://tree.opentreeoflife.org).

7. **Paleoclimate Data**: Paleoclimate data is not provided in this repository. It was obtained from the repository by Emily J. Judd et al. (https://github.com/EJJudd/PhanDA).

## Instructions for Use

To reproduce the analysis or apply these methods to your own dataset (e.g., a custom matrix of trait proportions per ecoregion), please follow these steps:

### 1. Data Preparation

Ensure your input data is formatted correctly.

* **For Main Analyses (Code01 - Code09, Code 19-23):** The scripts expect data aggregated at the ecoregion level (e.g., rows represent Ecoregions, columns represent trait proportions or environmental variables).
* **Preprocessing:** If you start with species-level data, use **Supplementary Code 1** to aggregate species occurrences into ecoregion trait proportions and perform Logit transformations.

### 2. Configuration (Path & Input Updates)

Each script contains a **Configuration** or **Parameter Settings** section at the beginning. You must update the working directory and input filenames to match your local environment.

**Example: Updating `Code01-PCA and Visualization of Logit-transformed Ecoregion Trait Proportions (Fig 1a, 1b).R`**

Open the script and locate Section 1 and 2. Modify the paths to point to your specific dataset:

```r
# 1. Set Working Directory
#-----------------------------------------------------------------------------
# Change this path to the folder containing your data
setwd("C:/Your/Local/Path/To/Data") 

# ... libraries loading ...

# 2. Parameter Settings
#-----------------------------------------------------------------------------
# --- File Paths ---
# REPLACE with your own aggregated ecoregion trait data file
input_data_path  <- "Your_Custom_Ecoregion_Trait_Data.csv" 

# Define your target output filenames
output_pca_data  <- "Output_PCA_Results.csv"
output_pdf_a     <- "Analysis_Fig1_a.pdf"
output_pdf_b     <- "Analysis_Fig1_b.pdf"
```

### 3. Running the Analysis

Run the individual scripts for specific analyses. You can execute them line-by-line in RStudio or source the entire file:

```r
# Example: Run PCA visualization
source("Code01-PCA and Visualization of Logit-transformed Ecoregion Trait Proportions (Fig 1a, 1b).R")

# Example: Run Spatial Piecewise SEM
source("Code23-Global Spatial Piecewise SEM (Fig 5).R")
```

### 4. Outputs

These scripts are designed to automatically save outputs (figures as PDF and processed data as CSV) to your working directory.

* Check the console for progress messages (e.g., `>>> Plot saved to...`).

## Reproduction Instructions
To fully reproduce the manuscript's results starting from raw files, please obtain the primary datasets from the following official sources:

 **Taxonomy**: The species list was derived from the **World Flora Online (WFO)** backbone taxonomy (Version 2025.02), available at [worldfloraonline.org](https://worldfloraonline.org).

 **Species Traits**: Primary trait data were extracted from digitized floras including eFloras(http://www.efloras.org), Flora of Australia(https://profiles.ala.org.au/opus/foa), and WFO. Cross-validation and supplementation were performed using TRY(https://www.try-db.org), GIFT(https://gift.uni-goettingen.de), and Tree of Sex(https://treeofsex.org).
 
 **Species Occurrences**: Occurrence records were obtained from the **Global Biodiversity Information Facility (GBIF)**. The specific dataset used in this study is available at **(https://doi.org/10.15468/dl.dj8zxb)**.
 
 **Ecoregion Boundaries**: Terrestrial ecoregions were delineated according to the **WWF** classification.
 
 **Environmental Data**:

   **Contemporary Climate (1979–2013)**: High-resolution raster data from **CHELSA V1.2** (http://www.chelsa-climate.org).
   **Paleoclimate (LGM, ~21 ka)**: Anomalies calculated using data from **Paleoclim V1.2b** (http://www.paleoclim.org).
   **Aridity & PET**: **Global Aridity Index and Potential Evapotranspiration Database v3** (https://www.global-ai-pet.org).
   **Solar Radiation, Wind, Elevation**: **WorldClim v2.1** (https://www.worldclim.org/).
   **Soil Properties**: Topsoil data (SOC, pH, Sand) from the **Harmonized World Soil Database v2.0 (HWSD)** (https://gaez.fao.org/pages/hwsd).

 **Phylogenetic Data**:

   Species-level modeling: **Smith and Brown backbone (GBOTB)** via Open Tree of Life(https://tree.opentreeoflife.org).
   Spatial-level analysis: `GBOTB.extended.TPL.tre` mega-tree provided in the **V.PhyloMaker2** R package.

Pre-processed datasets (e.g., aggregated ecoregion trait proportions, environmental data and principal components, and phylogenetic similarity matrices) are provided in the `Data` folder of this repository.

 **Run the scripts in sequence to reproduce**:

* `Supplementary Code1-Calculation of Ecoregion Trait Proportions and Logit-Transform.R`

* **Result Part 1: A Unified Global Spectrum of the Spatial Variation in Taxonomic Traits**
* `Code01-PCA and Visualization of Logit-transformed Ecoregion Trait Proportions (Fig 1a, 1b).R`

* **Result Part 2: Environmental drivers of the spatial variation spectrum**
* `Code02-PGLMM Analysis of Logit-Transformed Ecoregion Proportion.R`
* `Code03-PGLMM Analysis of PCA Axes (PC1-PC3).R`
* `Code04-Visualization of PGLMM Model Results (Fig 2a, 2b, 2c).R`
* `Code05-Data Process of Environmental PCA and Biplot Visualization (Fig 2d).R`
* `Code06-Mapping Trait PCs onto Environmental Space (GAM Smoothing).R`

* **Result Part 3: Phylogenetic imprints on spatial variation and species-level evolution**
* `Code07-The Second PGLMM Analysis of PCA Axes (PC1-PC3).R`
* `Code08-Visualization of PGLMM Coefficient (Fig 3a).R`
* `Code09-Phylo Metrics vs Trait PCs (Fig 3b).R`
* `Code10-Phylogenetic Signal of Binary Traits (Fritz's D via Subsampling).R`
* `Code11-Phylogenetic Conservation via consenTRAIT (Mean Trait Depth).R`
* `Code12-Evolutionary Model Fitting (ER vs ARD).R`
* `Code13-Summary Best Evolutionary Models for 15 Traits.R`
* `Code14-Stochastic Character Mappings Simulation.R`
* `Code15-Summarize SCM Simulation Results.R`
* `Code16-Beta Regression Analysis of Trait Proportions vs. Paleoclimate.R`
* `Code17-Data Integration of Trait Evolution Metrics.R`
* `Code18-Visualization Evolutionary Metrics (Fig 3c, d).R`

* **Result Part 4: Global patterns of taxonomic trait divergence**
* `Code19-Taxonomic Traits Diversity and Null Models.R`
* `Code20-Trait Divergence Mapping on Environmental Space (Fig 4b).R`
* `Code21-Trait Divergence Mapping on Trait Space (Fig 4d).R`
* `Code22-Trait Divergence and Trait Space (Fig 4c).R`

* **Result Part 5: Linking environment, phylogeny, and trait patterns**
* `Code23-Global Spatial Piecewise SEM (Fig 5).R`

**Expected run time**: Less than 1 minute per script for visualizations on a standard desktop computer; overall reproduction may take 4-24 hours depending on manual data preparation scale and hardware (parallelize loops for speed; e.g., use multicore for null models). The most time-intensive steps are the evolutionary modeling(Code12), the stochastic character mapping simulations (Code14) and null model simulations (Code19), which involve thousands of iterations. Please be aware that fitting these models may take several days to complete.

For questions, contact the code author (first author of the paper): Chunyu Zhang (chunyuzhanggeobio@gmail.com) or the corresponding authors: Hui Wang & Peiming Zheng (zhengpeiming@email.sdu.edu.cn).

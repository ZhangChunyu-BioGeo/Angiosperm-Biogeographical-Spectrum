#-----------------------------------------------------------------------------
# Code 9: Phylo Metrics vs Trait PCs (Fig 3b)
#
# Description:
# Generates a 2x4 matrix of hexbin plots showing the relationship between:
#   Rows: Trait PC1, PC2
#   Cols: Phylogenetic Metrics (PD, MDT, MPD, MNTD)
#
# Inputs:
#   - Global_Ecoregion_Phylogenetic.csv
#   - PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv
#
# Outputs:
#   - Fig3_b.pdf
#-----------------------------------------------------------------------------

# 0. Global Settings
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("data.table", "ggplot2", "hexbin", "cowplot", "grid")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(ggplot2)
library(hexbin)
library(cowplot)
library(grid)

# 1. Data Loading and Preparation
#-----------------------------------------------------------------------------
cat("1. Loading Data...\n")

phylo_path <- "Global_Ecoregion_Phylogenetic.csv"
trait_path <- "PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv"

if(!file.exists(phylo_path) || !file.exists(trait_path)) {
  stop("Error: Input files not found.")
}

# Load using data.table
df_phylo <- fread(phylo_path, select = c("ECO_ID", "PD", "MDT", "MPD", "MNTD"))
df_pca   <- fread(trait_path, select = c("ECO_ID", "PC1", "PC2"))

# Merge
df_merged <- merge(df_phylo, df_pca, by = "ECO_ID")
df_merged <- na.omit(df_merged)

cat(sprintf("   -> Analysis based on %d ecoregions.\n", nrow(df_merged)))

# 2. Visualization Settings
#-----------------------------------------------------------------------------
# Color Palettes
palette_pos <- c("#FDFBE9", "#D27EAC", "#8C1C54") # Cream -> Pink -> Purple (Positive)
palette_neg <- c("#FDFBE9", "#9DC86E", "#427638") # Cream -> Light Green -> Dark Green (Negative)
palette_ns  <- c("#E8E8E8", "#A6A6A6", "#7F7F7F") # Grey (Non-significant)

color_values <- c(0, 0.3, 1)

# Variables to Plot
row_vars <- c("PC1", "PC2")
col_vars <- c("PD", "MDT", "MPD", "MNTD")

# Layout Control
target_aspect_ratio <- 0.6  # Narrow width (Height = ~1.66 * Width)
plot_list <- list()

cat("2. Generating Hexbin Plots...\n")

# 3. Plotting Loop
#-----------------------------------------------------------------------------
for (r in seq_along(row_vars)) {
  for (c in seq_along(col_vars)) {
    
    y_name <- row_vars[r]
    x_name <- col_vars[c]
    
    # Extract Data for current pair
    sub_data <- df_merged[, .(X = get(x_name), Y = get(y_name))]
    
    # Calculate Statistics
    cor_res <- cor.test(sub_data$X, sub_data$Y)
    r_val   <- cor_res$estimate
    p_val   <- cor_res$p.value
    
    star <- dplyr::case_when(
      p_val < 0.001 ~ "***", 
      p_val < 0.01 ~ "**", 
      p_val < 0.05 ~ "*", 
      TRUE ~ "ns"
    )
    label_text <- paste0("R = ", sprintf("%.2f", r_val), " ", star)
    
    # Determine Colors
    if (p_val >= 0.05) {
      current_colors <- palette_ns
    } else if (r_val > 0) {
      current_colors <- palette_pos
    } else {
      current_colors <- palette_neg
    }
    
    # Labels Logic
    y_lab <- if (c == 1) y_name else NULL
    x_lab <- if (r == 2) x_name else NULL
    top_title <- if (r == 1) x_name else NULL
    
    # == Aspect Ratio Logic (Preserved from snippet) ==
    
    # 1. Data Range
    x_rng <- range(sub_data$X)
    y_rng <- range(sub_data$Y)
    x_span <- diff(x_rng)
    y_span <- diff(y_rng)
    
    # 2. Base Ratio for Regular Hexagons
    base_ratio <- x_span / y_span
    
    # 3. Expand Y-axis to achieve target narrowness
    required_y_span <- x_span / (base_ratio * target_aspect_ratio)
    
    # 4. New Y Limits
    y_mid <- mean(y_rng)
    y_half <- required_y_span / 2
    new_ylim <- c(y_mid - y_half, y_mid + y_half)
    
    # == Generate ggplot ==
    p <- ggplot(sub_data, aes(x = X, y = Y)) +
      
      # Hexbins (bins=10 for large hexes)
      geom_hex(bins = 10, color = "white", linewidth = 0.05) +
      
      # Trend Line
      geom_smooth(method = "lm", color = "black", fill = "gray85", 
                  alpha = 0.5, linewidth = 0.5, se = TRUE) +
      
      # Stats Label
      annotate("text", x = -Inf, y = Inf, label = label_text,
               hjust = -0.1, vjust = 1.3, 
               size = 3, color = "black", fontface = "plain") +
      
      # Lock Aspect Ratio
      coord_fixed(ratio = base_ratio, ylim = new_ylim) +
      
      scale_fill_gradientn(colors = current_colors, values = color_values, name = "Count") +
      labs(x = x_lab, y = y_lab, title = top_title) +
      
      theme_minimal(base_size = 11) +
      theme(
        text = element_text(face = "plain", color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank(),
        
        # Conditional Axes
        axis.text.y = if(c == 1) element_text(face = "plain", color = "black", size = 8) else element_blank(),
        axis.ticks.y = if(c == 1) element_line(color = "black", linewidth = 0.3) else element_blank(),
        axis.title.y = if(c == 1) element_text(face = "plain", size = 11, color = "black") else element_blank(),
        
        axis.text.x = if(r == 2) element_text(face = "plain", color = "black", size = 8) else element_blank(),
        axis.ticks.x = if(r == 2) element_line(color = "black", linewidth = 0.3) else element_blank(),
        axis.title.x = if(r == 2) element_text(face = "plain", size = 11, color = "black") else element_blank(),
        
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5, face = "plain", size = 11, margin = margin(b = 5)),
        legend.position = "none",
        
        plot.margin = margin(2, 2, 2, 2)
      )
    
    plot_list[[length(plot_list) + 1]] <- p
  }
}

# 4. Assemble and Save
#-----------------------------------------------------------------------------
cat("3. Saving Figure...\n")

final_plot <- plot_grid(
  plotlist = plot_list, 
  ncol = 4, 
  nrow = 2, 
  align = "hv", 
  axis = "lbtr"
)

print(final_plot)

# Save as Fig3_b.pdf (Following Fig3_a from Code 8)
ggsave("Fig3_b.pdf", final_plot, width = 12, height = 7)
cat(">>> Saved: Fig3_b.pdf\n")

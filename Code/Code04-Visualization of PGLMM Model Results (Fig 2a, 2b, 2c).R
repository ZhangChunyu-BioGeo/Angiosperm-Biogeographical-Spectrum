#-----------------------------------------------------------------------------
# Code 4: Visualization of PGLMM Model Results (Fig 2a, 2b, 2c)
#
# Description:
# 1. Fig 2a: Heatmap of PGLMM coefficients for 15 Traits.
# 2. Fig 2b: Forest plot of PGLMM coefficients for PC1-PC3.
# 3. Fig 2c: Variance partitioning stacked bar plot for PC1-PC3.
#
# Note: "SpeciesNumber" is visualized as "SpN".
#
# Inputs:
# - PGLMM_Coefficients_Standardized_Logit_Proportion.csv
# - PGLMM_Coefficients_Standardized_PCs.csv
# - PGLMM_Model_Statistics_PCs.csv
#
# Outputs:
# - Fig2_a.pdf
# - Fig2_b.pdf
# - Fig2_c.pdf
#-----------------------------------------------------------------------------

# 0. Global Settings and Library Loading
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("data.table", "pheatmap", "ggplot2", "grid", "gridExtra", "cowplot", "reshape2", "tidyr")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(pheatmap)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot) 
library(reshape2)
library(tidyr)

#=============================================================================
# PART 1: Fig 2a - PGLMM Heatmap with Standard ggplot Legend
#=============================================================================
cat("\n=== Preparing Fig 2a: PGLMM Heatmap ===\n")

# 1. Data Preparation
coef_file_path <- "PGLMM_Coefficients_Standardized_Logit_Proportion.csv"
if(!file.exists(coef_file_path)) stop("Error: Coefficients file not found!")
df_coef <- fread(coef_file_path)

# Filter target variables
target_env_vars <- c("AMT", "MDTR", "AP", "PS",  "SRAD", 
                     "Wind", "AMTd", "APd", "Elev", "Slope", "Sand", "SOC", "PH", 
                     "SpeciesNumber", "Area") # Input file uses SpeciesNumber
df_sub <- df_coef[Variable %in% target_env_vars]

# Rename SpeciesNumber -> SpN
df_sub[Variable == "SpeciesNumber", Variable := "SpN"]

# Construct Matrices
mat_beta <- dcast(df_sub, Variable ~ Trait, value.var = "Std_Beta")
rownames(mat_beta) <- mat_beta$Variable
mat_beta <- as.matrix(mat_beta[, -1])

mat_pval <- dcast(df_sub, Variable ~ Trait, value.var = "P_Value")
rownames(mat_pval) <- mat_pval$Variable
mat_pval <- as.matrix(mat_pval[, -1])

# Extract Clustered Order for Environmental Variables (for Part 2)
# We only want to cluster the environmental vars, keeping SpN and Area separate later
env_vars_only <- c("AMT", "MDTR", "AP", "PS",  "SRAD", 
                   "Wind", "AMTd", "APd", "Elev", "Slope", "Sand", "SOC", "PH")
mat_beta_env <- mat_beta[rownames(mat_beta) %in% env_vars_only, ]
dist_rows <- dist(mat_beta_env, method = "euclidean")
hclust_rows <- hclust(dist_rows, method = "complete")
env_order_core <- rownames(mat_beta_env)[hclust_rows$order]

# 2. Define Colors
custom_colors <- c("#427638", "#9DC86E", "#FDFBE9", "#D27EAC", "#8C1C54")
color_palette_pheatmap <- colorRampPalette(custom_colors)(100)

max_abs_val <- max(abs(mat_beta), na.rm = TRUE)
breaks_seq  <- seq(-max_abs_val, max_abs_val, length.out = 101)

# 3. Generate Main Heatmap (pheatmap)
heatmap_obj <- pheatmap(
  mat_beta,
  color = color_palette_pheatmap,
  breaks = breaks_seq,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  display_numbers = FALSE,
  angle_col = 45,
  fontsize = 10,
  border_color = "white",
  legend = FALSE,
  silent = TRUE
)

# --- Add Significance Circles ---
row_ord <- heatmap_obj$tree_row$order
col_ord <- heatmap_obj$tree_col$order
p_reordered <- mat_pval[row_ord, col_ord]

nr <- nrow(p_reordered); nc <- ncol(p_reordered)
x_coords <- matrix(rep((1:nc - 0.5) / nc, each = nr), nrow = nr)
y_coords <- matrix(rep((nr:1 - 0.5) / nr, times = nc), nrow = nr)

radii_matrix <- apply(p_reordered, c(1, 2), function(p) {
  cell_size_scaler <- min(1/nr, 1/nc) * 1.2 
  if (is.na(p) || p >= 0.05) { return(0) } 
  else if (p < 0.001) { return(cell_size_scaler * 0.3) } 
  else if (p < 0.01) { return(cell_size_scaler * 0.18) } 
  else { return(cell_size_scaler * 0.08) }
})

sig_indices <- which(radii_matrix > 0)
if (length(sig_indices) > 0) {
  circles_grob <- circleGrob(
    x = x_coords[sig_indices], y = y_coords[sig_indices],
    r = radii_matrix[sig_indices],
    gp = gpar(fill = "transparent", col = "black", lwd = 0.8) 
  )
  matrix_grob_index <- which(heatmap_obj$gtable$layout$name == "matrix")
  heatmap_obj$gtable$grobs[[matrix_grob_index]] <- 
    addGrob(heatmap_obj$gtable$grobs[[matrix_grob_index]], circles_grob)
}

# 4. Generate Dummy Plot for Standard ggplot Legend
dummy_data <- data.frame(val = seq(-max_abs_val, max_abs_val, length.out = 100), x=1, y=1)

p_dummy <- ggplot(dummy_data, aes(x=x, y=y, fill=val)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = custom_colors,
    limits = c(-max_abs_val, max_abs_val),
    name = "Std. Beta"
  ) +
  theme_void() + 
  theme(
    legend.title = element_text(face = "bold", size = 10, hjust = 0.5),
    legend.text = element_text(size = 9)
  ) +
  guides(fill = guide_colorbar(
    barwidth = unit(0.5, "cm"),   
    barheight = unit(5, "cm"),    
    ticks.colour = "black",       
    frame.colour = "black",       
    title.position = "top",       
    label.position = "right"      
  ))

legend_grob <- get_legend(p_dummy)

# 5. Assemble Fig 2a
p_fig2a <- plot_grid(
  heatmap_obj$gtable, 
  legend_grob,
  ncol = 2,
  rel_widths = c(1, 0.2)
)


#=============================================================================
# PART 2: Fig 2b - PGLMM Coefficient Forest Plot for PCs
#=============================================================================
cat("\n=== Preparing Fig 2b: PC Coefficients Forest Plot ===\n")

# 1. Prepare Labels (R2 Statistics)
file_stats <- "PGLMM_Model_Statistics_PCs.csv"
if(!file.exists(file_stats)) stop("Error: Stats file not found!")

df_stats <- fread(file_stats)
df_stats[, New_Label := sprintf("%s\n(R2m=%.3f, R2c=%.3f)", 
                                Trait, R2_m, R2_c)]
facet_labels <- setNames(df_stats$New_Label, df_stats$Trait)

# 2. Prepare Plotting Data
file_pcs <- "PGLMM_Coefficients_Standardized_PCs.csv"
if(!file.exists(file_pcs)) stop("Error: PC coefficients file not found!")

df_pcs <- fread(file_pcs)

# Define variables
df_pcs[Variable == "SpeciesNumber", Variable := "SpN"]

all_vars <- c(env_vars_only, "SpN", "Area")
df_plot <- df_pcs[Variable %in% all_vars]

# Apply Y-axis order
# Logic: Heatmap Order (Top) -> SpN -> Area (Bottom)
# ggplot builds from bottom up, so we reverse the list: Area -> SpN -> Heatmap Order
final_order_top_to_bottom <- c(env_order_core, "SpN", "Area")
df_plot$Variable <- factor(df_plot$Variable, levels = rev(final_order_top_to_bottom))

# 3. Aux Columns
df_plot$P_Cat <- cut(df_plot$P_Value, 
                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", "ns"))

ns_color <- "grey50" 
df_plot$Line_Color <- ifelse(df_plot$P_Value < 0.05, "black", ns_color)

# 4. Generate Forest Plot
max_beta_pc <- max(abs(df_plot$Std_Beta), na.rm = TRUE) * 1.05 

p_fig2b <- ggplot(df_plot, aes(x = Std_Beta, y = Variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = Std_Beta - 1.96 * Std_Error, 
                     xmax = Std_Beta + 1.96 * Std_Error,
                     color = Line_Color),
                 height = 0, linewidth = 0.3) +    
  geom_point(aes(fill = Std_Beta, size = P_Cat, color = Line_Color),
             shape = 21, stroke = 0.4) +           
  facet_grid(. ~ Trait, scales = "free_x", labeller = labeller(Trait = facet_labels)) +              
  scale_fill_gradientn(
    colors = custom_colors,
    limits = c(-max_beta_pc, max_beta_pc),
    name = "Std. Beta"
  ) +
  scale_color_identity() +             
  scale_size_manual(
    values = c("***" = 5, "**" = 3.5, "*" = 2.5, "ns" = 1.5),
    name = "Significance"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  labs(x = "Standardized Beta Coefficient (95% CI)", y = NULL) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey95"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, color = "black", lineheight = 1.1),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 9),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold")
  )


#=============================================================================
# PART 3: Fig 2c - PGLMM Variance Partitioning
#=============================================================================
cat("\n=== Preparing Fig 2c: Variance Partitioning ===\n")

# 1. Prepare Data
# Using df_stats loaded in Part 2
df_stats[, Total_Var := Var_Fixed + Var_Phylo + Var_Space + Var_Resid]

df_vp <- df_stats[, .(
  Trait,
  Fixed = (Var_Fixed / Total_Var) * 100,
  Phylo = (Var_Phylo / Total_Var) * 100,
  Space = (Var_Space / Total_Var) * 100,
  Resid = (Var_Resid / Total_Var) * 100
)]

df_vp_long <- pivot_longer(df_vp, cols = c("Fixed", "Phylo", "Space", "Resid"),
                           names_to = "Component", values_to = "Percentage")

# Stacking Order: Fixed (Bottom) -> Phylo -> Space -> Resid
df_vp_long$Component <- factor(df_vp_long$Component, 
                               levels = c("Resid", "Space", "Phylo", "Fixed"))

# 2. Colors & Labels
vp_colors <- c(
  "Fixed" = "#99C361", 
  "Phylo" = "#DB75A1", 
  "Space" = "#91C2EF", 
  "Resid" = "grey85"
)

# 3. Generate Plot
p_fig2c <- ggplot(df_vp_long, aes(x = Trait, y = Percentage, fill = Component)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = vp_colors,
    name = "Variance Component",
    breaks = c("Resid", "Space", "Phylo", "Fixed"), 
    labels = c("Residual", "Spatial", "Phylogeny", "Environment") # Requested labels
  ) +
  labs(y = "Variance Explained (%)", x = NULL) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10)
  )


#=============================================================================
# FINAL OUTPUT: Print and Save All Plots
#=============================================================================
cat("\n=== Saving Files ===\n")

# Figure 2a
print(p_fig2a)
ggsave("Fig2_a.pdf", p_fig2a, width = 8, height = 7)
cat("  -> Saved: Fig2_a.pdf\n")

# Figure 2b
print(p_fig2b)
ggsave("Fig2_b.pdf", p_fig2b, width = 10, height = 8)
cat("  -> Saved: Fig2_b.pdf\n")

# Figure 2c
print(p_fig2c)
ggsave("Fig2_c.pdf", p_fig2c, width = 6, height = 5)
cat("  -> Saved: Fig2_c.pdf\n")

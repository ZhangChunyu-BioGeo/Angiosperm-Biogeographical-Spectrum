#-----------------------------------------------------------------------------
# Code 8: Visualization of PGLMM Coefficient(Fig 3a)
#
# Description:
# Visualizes the standardized beta coefficients from the PGLMM model.
# Includes 95% Confidence Intervals and significance levels.
#
# Variables included: 
#   - Environment: EnvPC1, EnvPC2
#   - Phylogeny: PD, MDT, MPD
#   - Control: Area
#
# Inputs:
#   - PGLMM_Coefficients_Standardized_PCs_Phylo_Metric.csv
#   - PGLMM_Model_Statistics_PCs_Phylo_Metric.csv
#
# Outputs:
#   - Fig3_a.pdf
#-----------------------------------------------------------------------------

# 0. Global Settings
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("data.table", "ggplot2", "scales")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(ggplot2)
library(scales)

# 1. Prepare Facet Labels (R2 Statistics)
#-----------------------------------------------------------------------------
file_stats <- "PGLMM_Model_Statistics_PCs_Phylo_Metric.csv"
if(!file.exists(file_stats)) stop(paste("Error: File not found -", file_stats))

df_stats <- fread(file_stats)

# Create label: "PC1\n(R2m=0.xxx, R2c=0.xxx)"
df_stats[, New_Label := sprintf("%s\n(R2m=%.3f, R2c=%.3f)", 
                                Trait, R2_m, R2_c)]

# Named vector for labeller
facet_labels <- setNames(df_stats$New_Label, df_stats$Trait)

# 2. Prepare Plotting Data (Coefficients)
#-----------------------------------------------------------------------------
file_coefs <- "PGLMM_Coefficients_Standardized_PCs_Phylo_Metric.csv"
if(!file.exists(file_coefs)) stop(paste("Error: File not found -", file_coefs))

df_plot <- fread(file_coefs)

# Define Target Variables and Display Order (Top to Bottom)
target_order <- c("EnvPC1", "EnvPC2", "PD", "MDT", "MPD", "Area")

# Filter and set Factor Levels
# Note: ggplot draws y-axis from bottom up, so we reverse the order
df_plot <- df_plot[Variable %in% target_order]
df_plot$Variable <- factor(df_plot$Variable, levels = rev(target_order))

# 3. Aesthetics Setup
#-----------------------------------------------------------------------------
# A. Significance Categories (for point size)
df_plot$P_Cat <- cut(df_plot$P_Value, 
                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", "ns"))

# B. Significance Color (Black vs Grey)
ns_color <- "grey50" 
df_plot$Line_Color <- ifelse(df_plot$P_Value < 0.05, "black", ns_color)

# C. Color Palette (Consistent with previous figures)
custom_colors <- c("#427638", "#9DC86E", "#FDFBE9", "#D27EAC", "#8C1C54")
max_beta <- max(abs(df_plot$Std_Beta), na.rm = TRUE) * 1.05 

# 4. Generate Forest Plot
#-----------------------------------------------------------------------------
p <- ggplot(df_plot, aes(x = Std_Beta, y = Variable)) +
  
  # A. Vertical Reference Line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  
  # B. Error Bars (95% CI)
  geom_errorbarh(aes(xmin = Std_Beta - 1.96 * Std_Error, 
                     xmax = Std_Beta + 1.96 * Std_Error,
                     color = Line_Color),
                 height = 0,           
                 linewidth = 0.4) +    
  
  # C. Points (Fill = Coefficient, Size = Significance)
  geom_point(aes(fill = Std_Beta, size = P_Cat, color = Line_Color),
             shape = 21,               
             stroke = 0.5) +           
  
  # D. Faceting
  facet_grid(. ~ Trait, scales = "free_x", labeller = labeller(Trait = facet_labels)) +              
  
  # E. Scales
  scale_fill_gradientn(
    colors = custom_colors,
    limits = c(-max_beta, max_beta),
    name = "Std. Beta"
  ) +
  scale_color_identity() +             
  scale_size_manual(
    values = c("***" = 5, "**" = 3.5, "*" = 2.5, "ns" = 1.5),
    name = "Significance"
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  
  # F. Labels & Theme
  labs(x = "Standardized Beta Coefficient (95% CI)", y = NULL) +
  
  theme_minimal() +
  theme(
    # Borders & Grid
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey95"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Strip (Facet Headers)
    strip.background = element_blank(),
    strip.text = element_text(size = 11, color = "black", face = "bold", lineheight = 1.1),
    
    # Axes
    axis.text.y = element_text(color = "black", size = 10, face = "bold"),
    axis.text.x = element_text(color = "black", size = 9),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm"),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold")
  )

# 5. Output
#-----------------------------------------------------------------------------
print(p)

ggsave("Fig3_a.pdf", p, width = 10, height = 6)
cat(">>> Saved: Fig3_a.pdf\n")

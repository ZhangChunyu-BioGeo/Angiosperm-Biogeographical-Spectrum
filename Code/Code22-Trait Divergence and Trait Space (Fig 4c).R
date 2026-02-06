#-----------------------------------------------------------------------------
# Code 22: Trait Divergence and Trait Space (Fig 4c)
#
# Description:
#   Visualizes the relationship between TaxMPD (SES) and Trait Space PCs (1 & 2).
#   - PC1 Relationship: Fitted with a Quadratic Model (Vertex Form).
#   - PC2 Relationship: Fitted with a Linear Model.
#   - Color Scale: Colored by EnvPC1 using Quantile Breaks
#
# Inputs:
#   1. Diversity_Taxonomic_Traits_SES_Summary.csv (SES Values)
#   2. PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv (Trait PCs)
#   3. Global_Ecoregion_Environment_EnvPCs.csv (EnvPC1 for coloring)
#
# Outputs:
#   - Fig4_c.pdf
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
rm(list = ls())
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

# [Color Palette]
CUSTOM_COLS <- c("#427638", "#9DC86E", "#FDFBE9", "#D27EAC", "#8C1C54")

# Load Packages
pkgs <- c("data.table", "ggplot2", "patchwork", "scales")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

cat(">>> 1. Loading and Merging Data...\n")

# 2. Data Loading
#-----------------------------------------------------------------------------
# (A) TaxMPD_SES
df_ses <- fread("Diversity_Taxonomic_Traits_SES_Summary.csv", select = c("ECO_ID", "TaxMPD_SES"))

# (B) Trait PCA (PC1, PC2)
df_pca <- fread("PCA_std_Trait_Logit_Proportion_Ecoregion_No_Imputation.csv", select = c("ECO_ID", "PC1", "PC2"))

# (C) EnvPC1
df_env <- fread("Global_Ecoregion_Environment_EnvPCs.csv", select = c("ECO_ID", "EnvPC1"))

# 3. Merge Data
#-----------------------------------------------------------------------------
data_list <- list(df_ses, df_pca, df_env)
plot_data <- Reduce(function(x, y) merge(x, y, by = "ECO_ID", all = FALSE), data_list)
plot_data <- na.omit(plot_data)

cat(sprintf("   -> Data loaded: %d ecoregions.\n", nrow(plot_data)))

# 4. Prepare Quantile Color Scale
#-----------------------------------------------------------------------------
# Calculate breaks based on quantiles to ensure even color distribution
probs <- seq(0, 1, length.out = 11) 
q_breaks <- unique(quantile(plot_data$EnvPC1, probs = probs))

cat("   -> Color scale breaks calculated based on EnvPC1 quantiles.\n")

# 5. Statistical Modeling & Label Generation
#-----------------------------------------------------------------------------
cat(">>> 2. Modeling Relationships...\n")

# --- Model 1: PC1 (Quadratic) ---
# Fit: y = c + bx + ax^2
m1 <- lm(TaxMPD_SES ~ poly(PC1, 2, raw = TRUE), data = plot_data)
sum_m1 <- summary(m1)

# Extract Coefficients
coefs_m1 <- coef(m1)
a_val <- coefs_m1[3] # x^2 coefficient
b_val <- coefs_m1[2] # x coefficient
c_val <- coefs_m1[1] # Intercept

# Calculate Vertex (Symmetry Axis)
# h = -b / 2a
axis_x <- -b_val / (2 * a_val)
# k = c - b^2 / 4a (Vertex Y)
vertex_y <- c_val - (b_val^2) / (4 * a_val)

# Format Formula: y = a(x - h)^2 + k
sign_h <- ifelse(axis_x >= 0, "-", "+") 
sign_k <- ifelse(vertex_y >= 0, "+", "-")

formula_pc1 <- sprintf("y = %.2f(x %s %.2f)² %s %.2f", 
                       a_val, sign_h, abs(axis_x), sign_k, abs(vertex_y))

# Stats
r2_pc1 <- sum_m1$adj.r.squared
p_val_pc1 <- pf(sum_m1$fstatistic[1], sum_m1$fstatistic[2], sum_m1$fstatistic[3], lower.tail = FALSE)
p_str_pc1 <- ifelse(p_val_pc1 < 0.001, "P < 0.001", sprintf("P = %.3f", p_val_pc1))

label_pc1 <- sprintf("%s\nAdj. R² = %.3f\n%s", formula_pc1, r2_pc1, p_str_pc1)
cat(sprintf("   -> PC1 Model: %s\n", formula_pc1))


# --- Model 2: PC2 (Linear) ---
# Fit: y = mx + b
m2 <- lm(TaxMPD_SES ~ PC2, data = plot_data)
sum_m2 <- summary(m2)
coefs_m2 <- coef(m2)
slope <- coefs_m2[2]
intercept <- coefs_m2[1]

# Format Formula
sign_int <- ifelse(intercept >= 0, "+", "-")
formula_pc2 <- sprintf("y = %.2fx %s %.2f", slope, sign_int, abs(intercept))

# Stats
r2_pc2 <- sum_m2$adj.r.squared
p_val_pc2 <- pf(sum_m2$fstatistic[1], sum_m2$fstatistic[2], sum_m2$fstatistic[3], lower.tail = FALSE)
p_str_pc2 <- ifelse(p_val_pc2 < 0.001, "P < 0.001", sprintf("P = %.3f", p_val_pc2))

label_pc2 <- sprintf("%s\nAdj. R² = %.3f\n%s", formula_pc2, r2_pc2, p_str_pc2)
cat(sprintf("   -> PC2 Model: %s\n", formula_pc2))

# 6. Define Theme
#-----------------------------------------------------------------------------
clean_theme <- theme_minimal() +
  theme(
    # Lines
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    
    # Text
    axis.title = element_text(size = 12, face = "plain"),
    axis.text = element_text(size = 10, color = "black"),
    
    # Grid & Background
    panel.grid = element_blank(),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )

# 7. Plotting
#-----------------------------------------------------------------------------
cat(">>> 3. Generating Plots...\n")

# --- Plot 1: PC1 (Quadratic) ---
p1 <- ggplot(plot_data, aes(x = PC1, y = TaxMPD_SES)) +
  # Symmetry Axis Line
  geom_vline(xintercept = axis_x, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  
  # Points
  geom_point(aes(fill = EnvPC1), 
             shape = 21,          
             color = "white",     
             stroke = 0.2,        
             size = 2.5, 
             alpha = 0.9) +
  
  # Trend Line
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), 
              color = "black", linewidth = 0.5, 
              fill = "grey20", alpha = 0.15) +
  
  # Statistics Label
  annotate("text", x = min(plot_data$PC1), y = max(plot_data$TaxMPD_SES), 
           label = label_pc1, hjust = 0, vjust = 1, size = 3.5) +
  
  # Color Scale (Quantile)
  scale_fill_gradientn(
    colors = CUSTOM_COLS,
    values = scales::rescale(q_breaks), 
    breaks = q_breaks,                  
    labels = function(x) sprintf("%.1f", x),
    name = "EnvPC1",
    guide = guide_colorsteps(           
      barwidth = unit(0.4, "cm"),
      barheight = unit(5, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.3,
      ticks = TRUE,
      show.limits = TRUE
    )
  ) +
  labs(x = "PC1 (Trait)", y = "TaxMPD (SES)") +
  clean_theme

# --- Plot 2: PC2 (Linear) ---
p2 <- ggplot(plot_data, aes(x = PC2, y = TaxMPD_SES)) +
  # Points
  geom_point(aes(fill = EnvPC1), 
             shape = 21, 
             color = "white", 
             stroke = 0.2, 
             size = 2.5, 
             alpha = 0.9) +
  
  # Trend Line
  geom_smooth(method = "lm", formula = y ~ x, 
              color = "black", linewidth = 0.5, 
              fill = "grey20", alpha = 0.15) +
  
  # Statistics Label
  annotate("text", x = min(plot_data$PC2), y = max(plot_data$TaxMPD_SES), 
           label = label_pc2, hjust = 0, vjust = 1, size = 3.5) +
  
  # Color Scale (Shared)
  scale_fill_gradientn(
    colors = CUSTOM_COLS,
    values = scales::rescale(q_breaks),
    breaks = q_breaks,
    labels = function(x) sprintf("%.1f", x),
    name = "EnvPC1",
    guide = guide_colorsteps(
      barwidth = unit(0.4, "cm"),
      barheight = unit(5, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.3,
      ticks = TRUE,
      show.limits = TRUE
    )
  ) +
  labs(x = "PC2 (Trait)", y = NULL) +
  clean_theme

# 8. Combine and Save
#-----------------------------------------------------------------------------
# Collect guides ensures one common legend
final_plot <- p1 + p2 + plot_layout(guides = "collect")

output_file <- "Fig4_c.pdf"

# Save as PDF
ggsave(output_file, final_plot, width = 10, height = 5, device = cairo_pdf)

cat("========================================================\n")
cat(sprintf("Success! Figure saved to:\n%s\n", file.path(getwd(), output_file)))
cat("========================================================\n")

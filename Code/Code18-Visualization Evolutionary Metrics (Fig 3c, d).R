#-----------------------------------------------------------------------------
# Code 18: Visualization Evolutionary Metrics (Fig 3c, d)
#
# Description:
#   Generates a composite figure with two panels:
#   1. Left (Fig 3c): Lollipop chart of GMST Beta Coefficients (Climate Response),
#      sorted by Phylogenetic Signal (Fritz's D).
#   2. Right (Fig 3d): Scatter plot of Fritz's D vs. Evolutionary Depth,
#      showing the relationship between signal strength and clade age.
#
#
# Inputs:
#   - Trait_Evolution_Summary.csv
#
# Outputs:
#   - Fig3_cd.pdf
#-----------------------------------------------------------------------------

# 1. Global Settings and Libraries
#-----------------------------------------------------------------------------
setwd("F:/GitHub/ZhangChunyu-GeoBio/Angiosperm-Biogeographical-Spectrum/Data")

pkgs <- c("ggplot2", "dplyr", "ggrepel", "gridExtra")
for(p in pkgs) { if(!require(p, character.only = TRUE)) install.packages(p) }

library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)

# 2. Data Preparation
#-----------------------------------------------------------------------------
cat("1. Loading Data...\n")
if(!file.exists("Trait_Evolution_Summary.csv")) stop("Summary file not found.")

df <- read.csv("Trait_Evolution_Summary.csv", stringsAsFactors = FALSE)

# Feature Engineering
plot_data <- df %>%
  # Calculate 95% Confidence Intervals for Beta Coefficients
  mutate(
    CI_Min = GMST_Coef - 1.96 * GMST_Coef_sd,
    CI_Max = GMST_Coef + 1.96 * GMST_Coef_sd
  ) %>%
  # Sorting Logic: 
  # 1. Group by Trait
  # 2. Calculate mean Fritz D per trait
  # 3. Reorder Trait factor based on this mean (for Y-axis sorting)
  group_by(Trait) %>%
  mutate(Sort_Key = mean(Fritz_D, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Trait = reorder(Trait, Sort_Key)) %>%
  # Ensure numeric consistency
  mutate(consenTRAIT_Depth = as.numeric(consenTRAIT_Depth))

# [Crucial Step] Calculate Global Limits for Size Scale
# This ensures that a rate of 0.5 looks the same size in both plots.
tr_limits <- range(plot_data$Transition_Rate, na.rm = TRUE)
common_size_range <- c(2, 8) # Bubble size range (min, max)

cat(sprintf("   - Global Transition Rate Range: %.3f to %.3f\n", tr_limits[1], tr_limits[2]))

# 3. Define Graphical Theme
#-----------------------------------------------------------------------------
clean_theme <- theme_minimal() +
  theme(
    # Typography
    text = element_text(family = "sans", face = "plain", color = "black", size = 11),
    axis.title = element_text(face = "plain", size = 10),
    axis.text = element_text(face = "plain", color = "black", size = 9),
    
    # Borders & Lines
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    
    # Legend Styling
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.8, "cm"), # Increased size for clear bubbles
    legend.key = element_rect(fill = NA, color = NA),
    
    # Background
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# 4. Plot 1: GMST Response (Lollipop Chart)
#-----------------------------------------------------------------------------
cat("2. Generating Plot 1 (GMST Response)...\n")

p1 <- ggplot(plot_data, aes(y = Trait, x = GMST_Coef)) +
  # A. Vertical Zero Line (Neutral Response)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", size = 0.3) +
  
  # B. Connector Lines (Lollipop stick)
  geom_line(aes(group = Trait), color = "gray85", size = 0.6) +
  
  # C. Error Bars (95% CI)
  geom_errorbarh(aes(xmin = CI_Min, xmax = CI_Max), 
                 height = 0, color = "black", size = 0.3) +
  
  # D. Data Points (Bubbles)
  # Fill = Coefficient Strength, Size = Evolutionary Rate
  geom_point(aes(fill = GMST_Coef, size = Transition_Rate), 
             shape = 21, color = "black", stroke = 0.2, alpha = 1) +
  
  # E. Labels (Specific State Names, e.g., "GF_H")
  geom_text_repel(aes(label = Trait_Type), 
                  size = 2.8, 
                  color = "black",
                  segment.color = NA,
                  box.padding = 0.2) +
  
  # F. Scales
  scale_fill_gradient2(
    low = "#3781B0", mid = "white", high = "#D27EAC", midpoint = 0,
    name = "GMST Coef"
  ) +
  scale_size_continuous(
    limits = tr_limits, 
    range = common_size_range, 
    name = "Evol. Rate"
  ) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_discrete(expand = expansion(mult = 0.05)) +
  
  # G. Labels & Theme
  labs(x = "GMST Beta Coefficient", y = "Trait (Sorted by Fritz D)") +
  clean_theme

# 5. Plot 2: Evolutionary Depth vs. Signal
#-----------------------------------------------------------------------------
cat("3. Generating Plot 2 (Depth vs Signal)...\n")

# Prep A: Layering (Plot larger points behind smaller points for visibility)
p2_layer_data <- plot_data %>%
  arrange(desc(Transition_Rate))

# Prep B: Labels & Error Bars
# Since Fritz D and Depth are Trait-level metrics, they are duplicated for state 1/2.
# We distinct() them to avoid plotting the same error bar twice.
p2_label_data <- plot_data %>%
  select(Trait, consenTRAIT_Depth, Fritz_D, Fritz_D_sd) %>%
  distinct()

p2 <- ggplot() +
  
  # A. Trend Line (Loess)
  geom_smooth(data = p2_layer_data, 
              aes(x = consenTRAIT_Depth, y = Fritz_D),
              method = "loess", color = "gray80", fill = "gray90", 
              alpha = 0.3, size = 0.5) +
  
  # B. Error Bars (Fritz D SD) - Plotted only once per trait
  geom_errorbar(data = p2_label_data,
                aes(x = consenTRAIT_Depth, 
                    ymin = Fritz_D - Fritz_D_sd, 
                    ymax = Fritz_D + Fritz_D_sd),
                width = 0,       
                color = "black", 
                size = 0.3,
                alpha = 1) +
  
  # C. Data Points
  # Fill = Fritz D (Consistency), Size = Transition Rate
  geom_point(data = p2_layer_data, 
             aes(x = consenTRAIT_Depth, y = Fritz_D, 
                 fill = Fritz_D, size = Transition_Rate), 
             shape = 21, color = "black", stroke = 0.2, alpha = 1) +
  
  # D. Scales
  scale_x_reverse(expand = expansion(mult = 0.05)) + # Newer (0) on right
  
  scale_fill_gradientn(
    colors = c("#427638", "#9DC86E", "#FDFBE9", "#D27EAC", "#8C1C54"),
    name = "Fritz's D"
  ) +
  
  scale_size_continuous(
    limits = tr_limits, 
    range = common_size_range, 
    name = "Evol. Rate",
    breaks = pretty(tr_limits, n = 4)
  ) +
  
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  
  # E. Labels (Trait Names)
  geom_text_repel(data = p2_label_data,
                  aes(x = consenTRAIT_Depth, y = Fritz_D, label = Trait), 
                  size = 3, 
                  color = "gray30",
                  segment.color = NA,
                  max.overlaps = 20) +
  
  # F. Labels & Theme
  labs(x = "Evolutionary Depth (Ma)", y = "Fritz's D (Phylogenetic Signal)") +
  clean_theme +
  guides(
    # Custom Legend: Make size legend gray to distinguish from color scale
    size = guide_legend(order = 2, override.aes = list(fill = "gray70", color = "black")),
    fill = guide_colorbar(order = 1)
  )

# 6. Save Output
#-----------------------------------------------------------------------------
cat("4. Saving Figure...\n")

output_pdf <- "Fig3_cd.pdf"

# Create a PDF device
pdf(output_pdf, width = 14, height = 7) # Wide format for side-by-side
  grid.arrange(p1, p2, ncol = 2, widths = c(1, 1.1))
dev.off()

cat("Figure saved to:\n")
cat(file.path(getwd(), output_pdf), "\n")

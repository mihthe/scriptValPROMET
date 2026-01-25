### fig.2B - stratified random subset for clear visualization
# run simulation file
# run analysis file

# libraries
library(tidyverse)
library(scales)

# differences
difference <- lmed - sim$m_sur
abs_difference <- abs(difference)

# comprehensive data frame
plot_data <- data.frame(
  patient_id = 1:nrow(sim),
  difference = difference,
  abs_difference = abs_difference,
  surface = sim$Su,
  m_sur = sim$m_sur,
  lmed = lmed,
  location = sim$L,
  tumor_size = sim$Si
)

# order by tumor size (ascending)
plot_data_ordered <- plot_data %>%
  arrange(tumor_size)

# create stratified random subset 
set.seed(2024)  # For reproducibility
n_display <- 100  # Number of cases to display

# stratified sampling to maintain location proportions
plot_data_subset <- plot_data_ordered %>%
  group_by(location) %>%
  slice_sample(prop = n_display / nrow(plot_data_ordered)) %>%
  ungroup() %>%
  arrange(tumor_size)  # Re-order by tumor size after sampling

# reassign x positions for the subset
plot_data_subset$x_position <- 1:nrow(plot_data_subset)

# # print sampling summary
# cat("=== SAMPLING SUMMARY ===\n")
# cat("Original N:", nrow(plot_data_ordered), "\n")
# cat("Subset N:", nrow(plot_data_subset), "\n\n")
# cat("Location distribution:\n")
# print(table(plot_data_subset$location))
# cat("\nProportions preserved:\n")
# print(round(prop.table(table(plot_data_subset$location)), 3))

# SETUP
# location labels
location_names <- c("Colon-rectum", "Duodenum", "Small intestine", "Stomach", "Esophageal")
plot_data_subset$location_name <- location_names[plot_data_subset$location]

# colors
location_colors <- c(
  "#FFB6C1",   # 1: Colon-rectum - light pink/salmon
  "#ADD8E6",   # 2: Duodenum - light blue
  "#90EE90",   # 3: Small intestine - light green
  "#DDA0DD",   # 4: Stomach - light purple/plum
  "#FFDAB9"    # 5: Esophageal - peach/light orange
)

# segment colors
over_col <- "#2ECC71"
under_col <- "#E74C3C"

# FIGURE
jpeg(paste0("output/figures/", "Difference_by_Size_subset.jpg"), 
     units = "in", 
     width = 18, height = 12, res = 400)

# layout: main plot (top) + legend panel (bottom)
layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1))

# PANEL 1: main plot
par(mar = c(6, 8, 4, 4))

plot(plot_data_subset$x_position, plot_data_subset$difference, 
     type = "n",
     xlab = "", 
     ylab = "Difference (Predicted λ - True Mitotic Count)",
     main = "Prediction Error Ordered by Tumor Size (Ascending)",
     bty = "n",
     ylim = range(plot_data_subset$difference) * 1.15,
     xaxt = "n",
     cex.lab = 1.5,
     cex.main = 1.7,
     cex.axis = 1.3)

# vertical background bars - now visible with fewer cases!
for (i in 1:nrow(plot_data_subset)) {
  rect(xleft = plot_data_subset$x_position[i] - 0.5,
       xright = plot_data_subset$x_position[i] + 0.5,
       ybottom = par("usr")[3],
       ytop = par("usr")[4],
       col = alpha(location_colors[plot_data_subset$location[i]], 0.5),  # Higher alpha
       border = NA)
}

# horizontal reference lines
abline(h = seq(floor(min(plot_data_subset$difference)/10)*10, 
               ceiling(max(plot_data_subset$difference)/10)*10, by = 10),
       col = alpha("gray50", 0.3), lty = 1, lwd = 0.5)

# zero line
abline(h = 0, lty = 1, col = "black", lwd = 1.5)

# segments - thicker now since fewer cases
segments(x0 = plot_data_subset$x_position, 
         y0 = 0, 
         y1 = plot_data_subset$difference,
         col = ifelse(plot_data_subset$difference > 0, over_col, under_col),
         lwd = 2.5)  # Thicker for visibility

# X-axis: tumor sizes
axis(1, at = plot_data_subset$x_position, 
     labels = round(plot_data_subset$tumor_size, 1),
     las = 2, 
     cex.axis = 0.5,  # Larger text now
     col.axis = "darkgreen", 
     col.ticks = "darkgreen",
     line = 0, 
     tck = -0.015)

mtext("Tumor Size (standardized)", side = 1, line = 4.5, cex = 1.2, col = "darkgreen")

# sample size annotation
mtext(paste0("Stratified random sample: n = ", nrow(plot_data_subset), 
             " of ", nrow(plot_data_ordered), " total cases"),
      side = 3, adj = 1, cex = 1.0, line = 0.5, font = 3)  # Italic

# PANEL 2: legends
par(mar = c(1, 8, 1, 4))

plot.new()

# location legend (sx) - horizontal layout
legend("left",
       title = expression(bold("Tumor Location")),
       legend = location_names,
       pch = 22,
       pt.bg = alpha(location_colors, 0.5),
       col = "black",
       pt.cex = 2.5,
       cex = 1.2,
       bty = "n",
       horiz = FALSE,
       ncol = 5)

# prediction direction legend (dx)
legend("right",
       title = expression(bold("Prediction Direction")),
       legend = c("Overestimation (λ > y)", 
                  "Underestimation (λ < y)", 
                  "Perfect prediction"),
       col = c(over_col, under_col, "black"),
       lwd = 4, 
       cex = 1.2, 
       bty = "n")

dev.off()

# reset layout
layout(1)
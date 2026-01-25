### fig.2B
## like the posterior predictive check plot
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
  patient_id = val[-217],
  difference = difference,
  abs_difference = abs_difference,
  surface = validation$Su[-217],  m_sur = validation$m_sur[-217],
  lmed = lmed[-217],
  location = validation$L[-217],
  tumor_size = validation$Si[-217]
)



# order ONLY by tumor size (ascending)
plot_data_ordered <- plot_data %>%
  arrange(tumor_size)

# sequential x-position
plot_data_ordered$x_position <- 1:nrow(plot_data_ordered)

# location labels
location_names <- c("Colon-rectum", "Duodenum", "Small intestine", "Stomach", "Esophageal")
plot_data_ordered$location_name <- location_names[plot_data_ordered$location]

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

# create figure
jpeg(paste0("output/figures/", "Difference_by_Size_styled.jpg"), 
     units = "in", 
     width = 18, height = 12, res = 400)

# create layout: main plot on top (larger), legend panel at bottom (smaller)
layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1))

# PANEL 1: main plot
par(mar = c(6, 8, 4, 4))  # bottom, left, top, right

plot(plot_data_ordered$x_position, plot_data_ordered$difference, 
     type = "n",
     xlab = "", 
     ylab = "Difference (Predicted λ - True Mitotic Count)",
     main = "Prediction Error Ordered by Tumor Size (Ascending)",
     bty = "n",
     ylim = range(plot_data_ordered$difference) * 1.15,
     xaxt = "n",
     cex.lab = 1.5,
     cex.main = 1.7,
     cex.axis = 1.3)

# draw vertical background bars
for (i in 1:nrow(plot_data_ordered)) {
  rect(xleft = plot_data_ordered$x_position[i] - 0.5,
       xright = plot_data_ordered$x_position[i] + 0.5,
       ybottom = par("usr")[3],
       ytop = par("usr")[4],
       col = alpha(location_colors[plot_data_ordered$location[i]], 0.35),
       border = NA)
}

# horizontal reference lines
abline(h = seq(floor(min(plot_data_ordered$difference)/10)*10, 
               ceiling(max(plot_data_ordered$difference)/10)*10, by = 10),
       col = alpha("gray50", 0.3), lty = 1, lwd = 0.5)

# zero line
abline(h = 0, lty = 1, col = "black", lwd = 1.5)

# segments
segments(x0 = plot_data_ordered$x_position, 
         y0 = 0, 
         y1 = plot_data_ordered$difference,
         col = ifelse(plot_data_ordered$difference > 0, over_col, under_col),
         lwd = 1.5)

# X-axis
axis(1, at = plot_data_ordered$x_position, 
     labels = round(plot_data_ordered$tumor_size, 1),
     las = 2, cex.axis = 0.35,
     col.axis = "darkgreen", col.ticks = "darkgreen",
     line = 0, tck = -0.01)

mtext("Tumor Size (standardized)", side = 1, line = 4.5, cex = 1.2, col = "darkgreen")

# PANEL 2: legends
par(mar = c(1, 8, 1, 4))  

plot.new()  # empty plot for legends

# location legend (sx)
legend("left",
       title = expression(bold("Tumor Location")),
       legend = location_names,
       pch = 22,
       pt.bg = alpha(location_colors, 0.35),
       col = "black",
       pt.cex = 2.5,
       cex = 1.2,
       bty = "n",
       horiz = FALSE,
       ncol = 5)  # Horizontal layout

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

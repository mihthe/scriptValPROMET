### better plot but improve rectangles!!
## like the posterior predictive check plot

### fig.2B
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

# order ONLY by tumor size (ascending)
plot_data_ordered <- plot_data %>%
  arrange(tumor_size)

# sequential x-position
plot_data_ordered$x_position <- 1:nrow(plot_data_ordered)

# location labels
location_names <- c("Colon-rectum", "Duodenum", "Small intestine", "Stomach", "Esophageal")
plot_data_ordered$location_name <- location_names[plot_data_ordered$location]

# colors - matching your goal image style
location_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
over_col <- "#2ECC71"
under_col <- "#E74C3C"

### FIGURE: plot styled like posterior predictive check ###
jpeg(paste0("output/figures/", "Difference_by_Size.jpg"), 
     units = "in", 
     width = 18, height = 10, res = 400)

par(mar = c(8, 8, 4, 4))

plot(plot_data_ordered$x_position, plot_data_ordered$difference, 
     type = "n",
     xlab = "", 
     ylab = "Difference (Predicted λ - True Mitotic Count)",
     main = "Prediction Error Ordered by Tumor Size (Ascending)",
     bty = "l",
     ylim = range(plot_data_ordered$difference) * 1.15,
     xaxt = "n",
     cex.lab = 1.5,
     cex.main = 1.7,
     cex.axis = 1.3)

# Draw vertical FULL HEIGHT background bars for each case (like in goal image)
for (i in 1:nrow(plot_data_ordered)) {
  rect(xleft = plot_data_ordered$x_position[i] - 0.5,
       xright = plot_data_ordered$x_position[i] + 0.5,
       ybottom = par("usr")[3],  # bottom of plot area
       ytop = par("usr")[4],      # top of plot area
       col = alpha(location_colors[plot_data_ordered$location[i]], 0.2),
       border = NA)
}

# horizontal grid
abline(h = seq(floor(min(plot_data_ordered$difference)), 
               ceiling(max(plot_data_ordered$difference)), by = 5),
       col = "gray85", lty = 1)

# zero line (thicker and more prominent)
abline(h = 0, lty = 1, col = "black", lwd = 2)

# Draw the difference segments with directional coloring
segments(x0 = plot_data_ordered$x_position, 
         y0 = 0, 
         y1 = plot_data_ordered$difference,
         col = ifelse(plot_data_ordered$difference > 0, 
                      over_col, 
                      under_col),
         lwd = 2)

# # points on top (optional - can remove if you want cleaner look)
# points(plot_data_ordered$x_position, plot_data_ordered$difference, 
#        pch = 21,
#        bg = alpha(location_colors[plot_data_ordered$location], 0.8),
#        col = "black",
#        cex = 1,
#        lwd = 0.5)

# X-axis: Tumor sizes
axis(1, at = plot_data_ordered$x_position, 
     labels = round(plot_data_ordered$tumor_size, 1),
     las = 2,
     cex.axis = 0.5,
     col.axis = "darkgreen",
     col.ticks = "darkgreen",
     line = 0)
mtext("Tumor Size (standardized)", side = 1, line = 3, cex = 1.2, col = "darkgreen")

# Legend - prediction direction
legend("topright", 
       title = expression(bold("Prediction Direction")),
       legend = c("Overestimation (λ > y)", 
                  "Underestimation (λ < y)", 
                  "Perfect prediction"),
       col = c(over_col, under_col, "black"),
       lwd = 3, cex = 1.2, bty = "n")

# Location legend
legend("topleft",
       title = expression(bold("Tumor Location")),
       legend = location_names,
       fill = alpha(location_colors, 0.2),
       border = "black",
       cex = 1.1,
       bty = "n")

dev.off()
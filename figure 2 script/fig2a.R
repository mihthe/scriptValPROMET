### fig.2A 
# plot only 100 ptz
# lambda accuracy by available surface

## line overlapping!!

# libraries
library(tidyverse)
library(scales)
sim_sur_l <- rpois(n = length(lambda), lambda = lambda)
sim_sur <- matrix(sim_sur_l, nrow = nrow(lambda), ncol = ncol(lambda))

# raw and absolute differences
difference <- sim_sur_l - validation$m_sur

#difference <- difference
abs_difference <- abs(difference)


# comprehensive data frame
plot_data <- data.frame(
  patient_id = val,
  difference = difference,
  abs_difference = abs_difference,
  surface = validation$Su,  m_sur = validation$m_sur,
  lmed = lmed,
  location = validation$L,
  tumor_size = validation$Si
)

# only the first 100 patients (rows)
N=length(val)-1
plot_data_N <- plot_data[1:N, ]

# order by surface for first analysis
plot_data_surface <- plot_data_N[order(plot_data_N$surface), ]
plot_data_surface$x_position <- 1:nrow(plot_data_surface)

# # order by surface for first analysis (all patients)
# plot_data_surface <- plot_data[order(plot_data$surface), ]
# plot_data_surface$x_position <- 1:nrow(plot_data_surface)

# colors for locations
# maybe not (?)
location_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

### FIGURE 1: Difference by Surface ###
jpeg(paste0("output/figures/", "Difference_by_Surface.jpg"), 
     units = "in", 
     width = 20, height = 12, res = 400)

par(mfrow = c(2, 1), mar = c(10, 8, 6, 4), xaxs = "r", yaxs = "r")

### PANEL A: Raw difference (signed) ###
plot(plot_data_surface$x_position, plot_data_surface$difference, 
     type = "h",
     lwd = 2,
     col = ifelse(plot_data_surface$difference > 0, 
                  alpha("#2ECC71", 0.6), 
                  alpha("#E74C3C", 0.6)),
     xlab = "", 
     ylab = "Signed Difference (sim_surg - m_surg)",
     main = "Signed Prediction Error (Ordered by Biopsy Surface)",
     bty = "l",
     ylim = range(plot_data_surface$difference) * 1.15,
     xaxt = "n",
     cex.lab = 1.4,
     cex.main = 1.6,
     cex.axis = 1.2)

abline(h = 0, lty = 1, col = "black", lwd = 2)

# Add surface as continuous line 
#par(new = TRUE)
plot(plot_data_surface$x_position, plot_data_surface$surface,
     type = "l",
     lwd = 3,
     col = alpha("navyblue", 0.7),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "",
     ylim = c(0, max(plot_data_surface$surface) * 1.1),
     bty = "n")

# points(plot_data_surface$x_position, plot_data_surface$difference, 
#        pch = 21,
#        bg = alpha(location_colors[plot_data_surface$location], 0.7),
#        col = "black",
#        cex = 1.5)

# X-axes
# axis(1, at = plot_data_surface$x_position, 
#      labels = plot_data_surface$patient_id,
#      las = 2, cex.axis = 0.6)
# mtext("Patient ID", side = 1, line = 3.5, cex = 1.2)

axis(1, at = plot_data_surface$x_position, 
     labels = round(plot_data_surface$surface, 1),
     las = 2, cex.axis = 0.6,
     col.axis = "navyblue", line = 5)
mtext("Biopsy Surface (HPF)", side = 1, line = 7.5, cex = 1.2, col = "navyblue")

# legend("topleft",
#        legend = c("Overestimation", "Underestimation"),
#        col = c("#2ECC71", "#E74C3C"),
#        lwd = 3, cex = 1.2, bty = "n", inset = c(0, 0.3))  # Move up by 30% of the plot area 
# # to avoid overlapping
dev.off()

# ### PANEL B: Absolute difference (magnitude only) ###
# plot(plot_data_surface$x_position, plot_data_surface$abs_difference, 
#      type = "h",
#      lwd = 2,
#      col = alpha("steelblue", 0.6),
#      xlab = "", 
#      ylab = "Absolute Difference |λ - y|",
#      main = "B) Magnitude of Prediction Error (Ordered by Biopsy Surface)",
#      bty = "l",
#      ylim = c(0, max(plot_data_surface$abs_difference) * 1.15),
#      xaxt = "n",
#      cex.lab = 1.4,
#      cex.main = 1.6,
#      cex.axis = 1.2)

# # points(plot_data_surface$x_position, plot_data_surface$abs_difference, 
# #        pch = 21,
# #        bg = alpha(location_colors[plot_data_surface$location], 0.7),
# #        col = "black",
# #        cex = 1.5)

# # Add mean absolute error line
# abline(h = mean(plot_data_surface$abs_difference), 
#        col = "red", lwd = 2, lty = 2)

# # Add surface as continuous line 
# par(new = TRUE)
# plot(plot_data_surface$x_position, plot_data_surface$surface,
#      type = "l",
#      lwd = 3,
#      col = alpha("navyblue", 0.7),
#      xaxt = "n", yaxt = "n",
#      xlab = "", ylab = "",
#      ylim = c(0, max(plot_data_surface$surface) * 1.1),
#      bty = "n")

# # X-axes
# # axis(1, at = plot_data_surface$x_position, 
# #      labels = plot_data_surface$patient_id,
# #      las = 2, cex.axis = 0.6)
# # mtext("Patient ID", side = 1, line = 3.5, cex = 1.2)

# axis(1, at = plot_data_surface$x_position, 
#      labels = round(plot_data_surface$surface, 1),
#      las = 2, cex.axis = 0.6,
#      col.axis = "navyblue", line = 5)
# mtext("Biopsy Surface (HPF)", side = 1, line = 7.5, cex = 1.2, col = "navyblue")

# # legend("topleft",
# #        legend = c(paste0("MAE = ", round(mean(plot_data_surface$abs_difference), 2)),
# #                   "Location: see color"),
# #        col = c("red", NA),
# #        lwd = c(2, NA),
# #        lty = c(2, NA),
# #        pch = c(NA, 21),
# #        pt.bg = c(NA, "gray"),
# #        cex = 1.2, bty = "n")

# dev.off()
jpeg(paste0("output/figures/", "Posterior_Difference.jpg"), 
     units = "in", 
     width = 5, height = 5, res = 400)
## using the posterior
lambda_diff <- lambda - validation$m_sur
density(lambda_diff)
dens(lambda_diff, xlim = c(-15,15),   lwd = 3, 
main = 'Signed Prediction Error on Validation set', 
xlab = "Posterior Probability of Signed Difference (sim_surg - m_surg)", 
show.zero = TRUE, adj=0.5)
legend("topleft", legend = "model underestimate", bty="n")
legend("topright", legend = "model overestimate", bty="n")
dev.off()


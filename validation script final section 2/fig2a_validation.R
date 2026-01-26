### fig.2A - ADAPTED FOR VALIDATION DATA
# plot prediction error by available biopsy surface

# run validation file
# get validation results from the saved RDS file:
#validation_results <- readRDS("validation_output1/results/validation_results.rds")

# libraries
library(tidyverse)
library(scales)

# extract the validation data from validation_results$data
val_data <- validation_results$data

# -----------------------------------------------------------------------------
# prepare data for plotting
# -----------------------------------------------------------------------------

# calculate raw and absolute differences
# predicted_median is the model's prediction (λ)
# m_surg is the observed mitotic count on surgery (y)
difference <- val_data$predicted_median - val_data$m_surg
abs_difference <- abs(difference)

# comprehensive data frame for plotting
plot_data <- data.frame(
  patient_id = 1:nrow(val_data),
  difference = difference,
  abs_difference = abs_difference,
  surface = val_data$biopsy_surface_hpf,  # Biopsy surface in HPF
  m_surg = val_data$m_surg,               # Observed surgical count
  lmed = val_data$predicted_median,       # Predicted median (λ)
  location = val_data$L,                  # Tumor location
  tumor_size = val_data$size_mm           # Tumor size in mm
)

# check how many patients we have
n_patients <- nrow(plot_data)
cat("\nTotal patients for plotting:", n_patients, "\n")

# use all patients (recommended for validation)
plot_data_subset <- plot_data

# -----------------------------------------------------------------------------
# order by biopsy surface for the analysis
# -----------------------------------------------------------------------------

plot_data_surface <- plot_data_subset[order(plot_data_subset$surface), ]
plot_data_surface$x_position <- 1:nrow(plot_data_surface)

# colors for locations (matching validation file)
# 1 = Colon-rectum, 2 = Duodenum, 3 = Small intestine, 4 = Stomach
location_colors <- c("coral", "gold", "forestgreen", "steelblue")
location_names <- c("Colon-rectum", "Duodenum", "Small intestine", "Stomach")

# =============================================================================
# create figure
# =============================================================================

# create output directory if it doesn't exist
#dir.create("validation_output1/figures", showWarnings = FALSE, recursive = TRUE)

jpeg(paste0("validation_output1/figures/", "fig2a_Difference_by_Surface_VALIDATION.jpg"), 
     units = "in", 
     width = 20, height = 12, res = 400)

par(mfrow = c(2, 1), mar = c(10, 8, 6, 4), xaxs = "r", yaxs = "r")

# -----------------------------------------------------------------------------
# PANEL A: Raw difference (signed) — shows over/underestimation
# -----------------------------------------------------------------------------

plot(plot_data_surface$x_position, plot_data_surface$difference, 
     type = "h",
     lwd = 2,
     col = ifelse(plot_data_surface$difference > 0, 
                  alpha("#2ECC71", 0.6),   # Green for overestimation
                  alpha("#E74C3C", 0.6)),  # Red for underestimation
     xlab = "", 
     ylab = expression("Signed Difference (" * lambda * " - y)"),
     main = "A) Signed Prediction Error (Ordered by Biopsy Surface) - External Validation",
     bty = "l",
     ylim = range(plot_data_surface$difference) * 1.15,
     xaxt = "n",
     cex.lab = 1.4,
     cex.main = 1.6,
     cex.axis = 1.2)

# zero line (perfect prediction)
abline(h = 0, lty = 1, col = "black", lwd = 2)

# add surface as continuous line (secondary visualization)
par(new = TRUE)
plot(plot_data_surface$x_position, plot_data_surface$surface,
     type = "l",
     lwd = 3,
     col = alpha("navyblue", 0.7),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "",
     ylim = c(0, max(plot_data_surface$surface) * 1.1),
     bty = "n")

# X-axis: biopsy surface values
axis(1, at = plot_data_surface$x_position, 
     labels = round(plot_data_surface$surface, 1),
     las = 2, cex.axis = 0.6,
     col.axis = "navyblue", line = 5)
mtext("Biopsy Surface (HPF)", side = 1, line = 7.5, cex = 1.2, col = "navyblue")

# legend for Panel A
legend("topleft",
       legend = c("Overestimation (λ > y)", 
                  "Underestimation (λ < y)", 
                  "Biopsy Surface"),
       col = c("#2ECC71", "#E74C3C", alpha("navyblue", 0.7)),
       lwd = c(3, 3, 3),
       lty = c(1, 1, 1),
       cex = 1.0, bty = "n")

# -----------------------------------------------------------------------------
# PANEL B: Absolute difference (magnitude only)
# -----------------------------------------------------------------------------

plot(plot_data_surface$x_position, plot_data_surface$abs_difference, 
     type = "h",
     lwd = 2,
     col = alpha("steelblue", 0.6),
     xlab = "", 
     ylab = expression("Absolute Difference |" * lambda * " - y|"),
     main = "B) Magnitude of Prediction Error (Ordered by Biopsy Surface) - External Validation",
     bty = "l",
     ylim = c(0, max(plot_data_surface$abs_difference) * 1.15),
     xaxt = "n",
     cex.lab = 1.4,
     cex.main = 1.6,
     cex.axis = 1.2)

# mean absolute error line
mae_value <- mean(plot_data_surface$abs_difference)
abline(h = mae_value, col = "red", lwd = 2, lty = 2)

# surface as continuous line
par(new = TRUE)
plot(plot_data_surface$x_position, plot_data_surface$surface,
     type = "l",
     lwd = 3,
     col = alpha("navyblue", 0.7),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "",
     ylim = c(0, max(plot_data_surface$surface) * 1.1),
     bty = "n")

# X-axis: biopsy surface values
axis(1, at = plot_data_surface$x_position, 
     labels = round(plot_data_surface$surface, 1),
     las = 2, cex.axis = 0.6,
     col.axis = "navyblue", line = 5)
mtext("Biopsy Surface (HPF)", side = 1, line = 7.5, cex = 1.2, col = "navyblue")

# legend for Panel B
legend("topleft",
       legend = c(paste0("MAE = ", round(mae_value, 2), " mitoses"),
                  "Biopsy Surface"),
       col = c("red", alpha("navyblue", 0.7)),
       lwd = c(2, 3),
       lty = c(2, 1),
       cex = 1.0, bty = "n")

dev.off()

# =============================================================================
# SUMMARY STATISTICS if needed
# =============================================================================

# cat("\n", paste(rep("═", 60), collapse = ""), "\n")
# cat("SUMMARY STATISTICS FOR FIG 2A\n")
# cat(paste(rep("═", 60), collapse = ""), "\n\n")
# 
# cat("Number of patients:", nrow(plot_data_surface), "\n")
# cat("Biopsy surface range:", round(min(plot_data_surface$surface), 1), "-", 
#     round(max(plot_data_surface$surface), 1), "HPF\n")
# cat("\nPrediction Error Statistics:\n")
# cat("  Mean Absolute Error (MAE):", round(mean(abs_difference), 2), "mitoses\n")
# cat("  Median Absolute Error:", round(median(abs_difference), 2), "mitoses\n")
# cat("  Mean Signed Error (Bias):", round(mean(difference), 2), "mitoses\n")
# cat("  SD of Error:", round(sd(difference), 2), "mitoses\n")
# 
# # Count over/under estimations
# n_over <- sum(difference > 0)
# n_under <- sum(difference < 0)
# n_exact <- sum(difference == 0)
# 
# cat("\nPrediction Direction:\n")
# cat("  Overestimations (λ > y):", n_over, "(", round(100*n_over/n_patients, 1), "%)\n")
# cat("  Underestimations (λ < y):", n_under, "(", round(100*n_under/n_patients, 1), "%)\n")
# cat("  Exact predictions (λ = y):", n_exact, "(", round(100*n_exact/n_patients, 1), "%)\n")
# 
# # Correlation between surface and absolute error
# cor_surface_error <- cor(plot_data_surface$surface, plot_data_surface$abs_difference)
# cat("\nCorrelation (Surface vs |Error|):", round(cor_surface_error, 3), "\n")
# cat("  Interpretation: ", 
#     ifelse(cor_surface_error < -0.3, "Larger surface → smaller errors (good!)",
#            ifelse(cor_surface_error > 0.3, "Larger surface → larger errors (unexpected)",
#                   "No strong relationship")), "\n")

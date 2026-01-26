## okish but weird violins

### fig1 validation
# run validation script and save results
#results <- readRDS("validation_output1/results/validation_results.rds")

# =============================================================================

# libraries
library(scales)
library(rethinking)  

# output directory if it doesn't exist
#dir.create("validation_output1/figures", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# extract data for plotting
# =============================================================================

# validation data with predictions
validation_data <- results$data

# posterior samples matrix: rows = posterior samples, cols = patients
# each column contains the posterior predictive distribution for one patient
lambda <- results$posterior_samples

# n of patients
N <- nrow(validation_data)
n_samples <- nrow(lambda)

cat(sprintf("Patients: %d\n", N))
cat(sprintf("Posterior samples per patient: %d\n", n_samples))

# =============================================================================
# prepare data for plotting
# =============================================================================

# median lambda for each patient (for ordering)
lambda_median <- apply(lambda, 2, median)

# order patients by ascending median predicted value 
idx <- order(lambda_median)

# y-axis limits based on observed surgical counts and predictions
y_min <- 0
y_max <- max(c(validation_data$m_surg, 
               apply(lambda, 2, quantile, probs = 0.975))) * 1.05

# cap y_max at a reasonable value if there are extreme outliers
y_max <- min(y_max, 80)  # adjust this cap as needed

cat(sprintf("Y-axis range: %.1f to %.1f\n", y_min, y_max))

# =============================================================================
# create plot
# =============================================================================

jpeg("validation_output1/figures/fig1_validation.jpg", 
     units = "in", 
     width = 14, height = 8, res = 300)

par(mai = c(0.9, 0.9, 0.9, 0.4))

# initialize empty plot
plot(NULL, 
     xlim = c(0.5, N + 0.5),  
     ylim = c(y_min, y_max),  
     xlab = 'Patient (ordered by predicted mitotic count)', 
     ylab = 'Mitotic count (per 5 mm²)',
     main = 'Figure 1: Posterior Predictive Check - External Validation',
     bty = 'l',
     xaxt = 'n')  

# x-axis with patient numbers
axis(1, at = 1:N, labels = 1:N, las = 2, cex.axis = 0.7)

# reference line at mitotic count = 5 (clinical threshold)
abline(h = 5, lty = 2, col = "gray40", lwd = 2)
text(N * 0.02, 6, "Risk threshold", adj = 0, cex = 0.7, col = "gray40")

# scaling factor for violin width
s_factor <- 0.4  

# loop through each patient (in ordered sequence)
for(i in 1:N) {
  patient_idx <- idx[i]  # the actual patient index in original data
  
  # posterior distribution for this patient
  patient_lambda <- lambda[, patient_idx]
  
  # density of the posterior distribution
  dens_obj <- density(patient_lambda, adjust = 1.2, from = 0)
  x_dens <- dens_obj$y  # density values (determines width of violin)
  y_dens <- dens_obj$x  # lambda values (determines height of violin)
  
  # normalise density for plotting (to avoid overlap between violins)
  x_dens_norm <- x_dens / max(x_dens) * s_factor
  
  # dx half of violin
  polygon(i + x_dens_norm, y_dens, 
          col = scales::alpha("skyblue", 0.5), 
          border = scales::alpha("steelblue", 0.7), 
          lwd = 0.3)
  
  # sx half of violin (mirror image)
  polygon(i - x_dens_norm, y_dens, 
          col = scales::alpha("skyblue", 0.5), 
          border = scales::alpha("steelblue", 0.7), 
          lwd = 0.3)
  
  # observed surgical mitotic count (m_surg) as blue cross
  points(i, validation_data$m_surg[patient_idx], 
         pch = 4, col = "blue", lwd = 2, cex = 1.2)
  
  # median prediction as a horizontal black line
  segments(i - 0.35, lambda_median[patient_idx], 
           i + 0.35, lambda_median[patient_idx], 
           lwd = 2, col = "black")
}

# legend
legend('topleft', 
       pch = c(4, NA, 15),
       lty = c(NA, 1, NA),
       lwd = c(2, 2, NA),
       col = c("blue", "black", scales::alpha("skyblue", 0.7)),
       pt.cex = c(1.2, NA, 2),
       legend = c("Observed (surgery)", 
                  "Predicted (median)",
                  "Posterior distribution"),
       bty = "n")

# subtitle
mtext("Violin shapes show full posterior predictive distribution for each patient", 
      side = 3, line = 0.3, cex = 0.8, adj = 0)

dev.off()

# =============================================================================
# SUMMARY STATISTICS if needed
# =============================================================================

# # coverage: how often does the observed value fall within the 89% HPDI?
# coverage <- numeric(N)
# for(i in 1:N) {
#   hpdi <- HPDI(lambda[, i], prob = 0.89)
#   coverage[i] <- (validation_data$m_surg[i] >= hpdi[1]) & 
#     (validation_data$m_surg[i] <= hpdi[2])
# }
# 
# cat(sprintf("\nModel Performance Summary:\n"))
# cat(sprintf("  - 89%% HPDI Coverage: %.1f%%\n", mean(coverage) * 100))
# cat(sprintf("  - MAE: %.2f mitoses\n", mean(abs(lambda_median - validation_data$m_surg))))
# cat(sprintf("  - Correlation (predicted vs observed): %.3f\n", 
#             cor(lambda_median, validation_data$m_surg)))

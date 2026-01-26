### RESIDUAL DIAGNOSTIC PLOTS FOR EXTERNAL VALIDATION

# run validation and save results

# libraries
library(rethinking)
library(scales)  

# output directory if it doesn't exist
#dir.create("validation_output1/figures", showWarnings = FALSE, recursive = TRUE)

# extract the components we need
validation_data <- results$data
posterior_samples <- results$posterior_samples  # matrix: [n_samples x n_patients]

# =============================================================================
# prepare data
# =============================================================================

# true mitotic counts on surgery (observed values)
Y_t <- validation_data$m_surg

# n of observations (patients)
N <- length(Y_t)

# n of posterior samples
l_s <- nrow(posterior_samples)

# posterior_samples matrix contains predicted counts (already Poisson draws)
# transpose to match the original code structure: [N x l_s]
Y_p <- t(posterior_samples)

# =============================================================================
# COMPUTE RESIDUALS
# =============================================================================

# residuals matrix: true - predicted (for each posterior sample)
resid <- matrix(nrow = N, ncol = l_s)

for(i in 1:N) {
  resid[i, ] <- Y_t[i] - Y_p[i, ]
}

# global HPDI of all residuals
resid_hpdi <- HPDI(c(resid), prob = 0.89)

# per-observation summaries
resid_median <- apply(resid, 1, median)
resid_hpdi_per_obs <- apply(resid, 1, HPDI, prob = 0.89)

cat(sprintf("Global 89%% HPDI of residuals: [%.2f, %.2f]\n", 
            resid_hpdi[1], resid_hpdi[2]))
cat(sprintf("Mean of median residuals: %.2f\n", mean(resid_median)))

# =============================================================================
# PLOT 1: RESIDUAL DISTRIBUTION
# =============================================================================

jpeg("validation_output1/figures/Residuals_distribution_validation.jpg", 
     units = "in", 
     width = 8, height = 6, res = 300)

# appropriate x-axis limits
x_range <- range(c(resid))
x_lim <- c(max(-50, x_range[1]), min(50, x_range[2]))

hist(resid, freq = FALSE, 
     main = "Posterior Distribution of Residuals\n(External Validation)", 
     xlab = "Residual (Observed - Predicted)",
     col = scales::alpha("steelblue", 0.5),
     border = "white",
     breaks = 50,
     xlim = x_lim)
abline(v = 0, col = "red", lwd = 2, lty = 2)
abline(v = resid_hpdi, col = "darkblue", lwd = 2, lty = 3)
legend("topright", 
       legend = c("Zero error", "89% HPDI"),
       col = c("red", "darkblue"),
       lwd = 2, lty = c(2, 3))

dev.off()

# =============================================================================
# PLOT 2: RESIDUALS VS FITTED 
# =============================================================================

# lambda (expected value) from posterior samples
# use the predicted_mean from validation_data
lambda_median <- validation_data$predicted_median

jpeg("validation_output1/figures/Residuals_vs_Fitted_validation.jpg", 
     units = "in", 
     width = 10, height = 8, res = 300)
par(mfrow = c(2, 2))

# Panel 1: Residuals vs Fitted Values
plot(lambda_median, resid_median,
     xlab = "Fitted values (predicted median)",
     ylab = "Median residual",
     main = "Residuals vs Fitted Values",
     pch = 16, col = scales::alpha(1, 0.5))
abline(h = 0, col = "red", lwd = 2, lty = 2)

# LOESS smooth if enough data points
if(N >= 10) {
  smooth_fit <- loess(resid_median ~ lambda_median)
  ord <- order(lambda_median)
  lines(lambda_median[ord], predict(smooth_fit)[ord], 
        col = "blue", lwd = 2)
}

# Panel 2: Scale-Location Plot (for heteroscedasticity)
abs_resid_median <- apply(abs(resid), 1, median)
plot(lambda_median, abs_resid_median,
     xlab = "Fitted values (predicted median)",
     ylab = "Median absolute residual",
     main = "Scale-Location Plot",
     pch = 16, col = scales::alpha(1, 0.5))

if(N >= 10) {
  smooth_fit2 <- loess(abs_resid_median ~ lambda_median)
  ord <- order(lambda_median)
  lines(lambda_median[ord], predict(smooth_fit2)[ord], 
        col = "blue", lwd = 2)
}

# Panel 3: Q-Q Plot of residuals
qqnorm(resid_median, 
       main = "Q-Q Plot of Median Residuals",
       pch = 16, col = scales::alpha(1, 0.5))
qqline(resid_median, col = "red", lwd = 2)

# Panel 4: Histogram with normal overlay
hist(resid_median, freq = FALSE,
     main = "Distribution of Median Residuals",
     xlab = "Median Residual",
     col = scales::alpha("steelblue", 0.5),
     border = "white",
     breaks = 15)
curve(dnorm(x, mean(resid_median), sd(resid_median)),
      add = TRUE, col = "red", lwd = 2)
abline(v = 0, col = "darkred", lwd = 2, lty = 2)

dev.off()

# =============================================================================
# PLOT 3: OBSERVED VS PREDICTED
# =============================================================================

jpeg("validation_output1/figures/Observed_vs_Predicted_validation.jpg", 
     units = "in", 
     width = 10, height = 10, res = 300)
par(mfrow = c(1, 1))

# prediction intervals
pred_median <- apply(Y_p, 1, median)
pred_hpdi <- apply(Y_p, 1, HPDI, prob = 0.89)

# location colors (matching validation_and_figures_correct.R)
loc_colors <- c("coral", "gold", "forestgreen", "steelblue")
loc_names <- c("Colon-rectum", "Duodenum", "Small intestine", "Stomach")

# point colors based on location
point_colors <- loc_colors[validation_data$L]

# axis limits
max_val <- max(c(Y_t, pred_median, pred_hpdi), na.rm = TRUE) + 2
max_val <- min(max_val, 100)  # cap at 100 for readability

plot(Y_t, pred_median,
     xlim = c(0, max_val),
     ylim = c(0, max_val),
     xlab = "Observed mitotic count (surgery)",
     ylab = "Predicted mitotic count (median)",
     main = "Observed vs Predicted Values\nwith 89% Prediction Intervals (External Validation)",
     pch = 16, col = scales::alpha(point_colors, 0.7),
     cex = 1.2)

# prediction intervals
segments(x0 = Y_t, 
         y0 = pmin(pred_hpdi[1, ], max_val), 
         y1 = pmin(pred_hpdi[2, ], max_val),
         col = scales::alpha("gray40", 0.3))

# prediction line
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

# legend for locations
legend("topleft", 
       legend = loc_names,
       pch = 16,
       col = scales::alpha(loc_colors, 0.7),
       title = "Location",
       pt.cex = 1.2)

# correlation info
r_squared <- cor(Y_t, pred_median)^2
legend("bottomright",
       legend = c(sprintf("R² = %.3f", r_squared),
                  sprintf("N = %d", N)),
       bty = "n")

dev.off()

# =============================================================================
# PLOT 4: RESIDUALS BY COVARIATES
# =============================================================================

# external validation data without
#   - Unknown confounder U (only exists in simulations)
#   - NAC response (resp_Rx) may not be available
# plot available covariates

jpeg("validation_output1/figures/Residuals_by_Covariates_validation.jpg", 
     units = "in", 
     width = 12, height = 10, res = 300)  # Increased height for margin space

# check which covariates are available
has_resp_Rx <- "resp_Rx" %in% names(validation_data) || 
  "n" %in% names(validation_data)

if(has_resp_Rx) {
  par(mfrow = c(2, 3), mar = c(6, 4, 3, 2))  
} else {
  par(mfrow = c(2, 2), mar = c(6, 4, 3, 2))  
}

# Panel 1: Residuals vs Tumor Size
plot(validation_data$size_mm, resid_median,
     xlab = "Tumor Size (mm)",
     ylab = "Median Residual",
     main = "Residuals vs Tumor Size",
     pch = 16, col = scales::alpha(point_colors, 0.6))
abline(h = 0, col = "red", lwd = 2, lty = 2)

# LOESS if enough points
if(N >= 10) {
  smooth_size <- loess(resid_median ~ validation_data$size_mm)
  ord <- order(validation_data$size_mm)
  lines(validation_data$size_mm[ord], predict(smooth_size)[ord], 
        col = "blue", lwd = 2)
}

# Panel 2: Residuals vs Biopsy Surface
plot(validation_data$biopsy_surface_hpf, resid_median,
     xlab = "Biopsy Surface (HPF)",
     ylab = "Median Residual",
     main = "Residuals vs Biopsy Surface",
     pch = 16, col = scales::alpha(point_colors, 0.6))
abline(h = 0, col = "red", lwd = 2, lty = 2)

if(N >= 10) {
  smooth_surf <- loess(resid_median ~ validation_data$biopsy_surface_hpf)
  ord <- order(validation_data$biopsy_surface_hpf)
  lines(validation_data$biopsy_surface_hpf[ord], predict(smooth_surf)[ord], 
        col = "blue", lwd = 2)
}

# Panel 3: Residuals vs Mitosis on Biopsy
plot(validation_data$m_bio, resid_median,
     xlab = "Mitotic Count on Biopsy",
     ylab = "Median Residual",
     main = "Residuals vs Biopsy Count",
     pch = 16, col = scales::alpha(point_colors, 0.6))
abline(h = 0, col = "red", lwd = 2, lty = 2)

if(N >= 10) {
  smooth_bio <- loess(resid_median ~ validation_data$m_bio)
  ord <- order(validation_data$m_bio)
  lines(validation_data$m_bio[ord], predict(smooth_bio)[ord], 
        col = "blue", lwd = 2)
}

# Panel 4: Residuals by Location (boxplot)
# shorter names for better fit
loc_names_short <- c("Colon-rect", "Duodenum", "Small int", "Stomach")
boxplot(resid_median ~ validation_data$L,
        names = loc_names_short,
        xlab = "Location",
        ylab = "Median Residual",
        main = "Residuals by Location",
        col = scales::alpha(loc_colors, 0.5),
        las = 2,
        cex.axis = 0.9)  # Slightly smaller axis text
abline(h = 0, col = "red", lwd = 2, lty = 2)

# # Panel 5: Residuals by NAC Response (if available)
# if(has_resp_Rx) {
#   # Determine the correct column name
#   if("resp_Rx" %in% names(validation_data)) {
#     nac_var <- validation_data$resp_Rx
#   } else if("n" %in% names(validation_data)) {
#     # n = 1 - resp_Rx, so Response = 0 in n
#     nac_var <- 1 - validation_data$n
#   }
#   
#   boxplot(resid_median ~ nac_var,
#           names = c("No Response/No NAC", "NAC Response"),
#           xlab = "NAC Response",
#           ylab = "Median Residual",
#           main = "Residuals by Treatment Response",
#           col = scales::alpha(c("coral", "steelblue"), 0.5))
#   abline(h = 0, col = "red", lwd = 2, lty = 2)
#   
#   # Panel 6: Note about missing confounder
#   plot.new()
#   text(0.5, 0.5, 
#        "Note: Unknown confounder (U)\nis not available in\nexternal validation data\n\n(U only exists in simulations)",
#        cex = 1.2, col = "gray40")
# }

dev.off()

# =============================================================================
# SUMMARY STATISTICS if needed
# =============================================================================

# cat("\n")
# cat("═══════════════════════════════════════════════════════════════════════\n")
# cat("RESIDUAL ANALYSIS SUMMARY (External Validation)\n")
# cat("═══════════════════════════════════════════════════════════════════════\n")
# cat(sprintf("Number of patients: %d\n", N))
# cat(sprintf("Number of posterior samples: %d\n", l_s))
# cat("\n")
# cat("RESIDUAL STATISTICS (Observed - Predicted):\n")
# cat(sprintf("  Mean residual:        %+.2f mitoses\n", mean(resid_median)))
# cat(sprintf("  Median residual:      %+.2f mitoses\n", median(resid_median)))
# cat(sprintf("  SD of residuals:       %.2f mitoses\n", sd(resid_median)))
# cat(sprintf("  89%% HPDI:             [%.2f, %.2f]\n", resid_hpdi[1], resid_hpdi[2]))
# cat("\n")
# cat(sprintf("  Proportion underestimated: %.1f%%\n", 100 * mean(resid_median > 0)))
# cat(sprintf("  Proportion overestimated:  %.1f%%\n", 100 * mean(resid_median < 0)))
# cat(sprintf("  Proportion exact:          %.1f%%\n", 100 * mean(resid_median == 0)))
# cat("\n")
# 
# # Coverage: how often does the 89% HPDI contain the true value?
# in_hpdi <- mapply(function(i) {
#   Y_t[i] >= resid_hpdi_per_obs[1, i] + Y_t[i] & 
#     Y_t[i] <= resid_hpdi_per_obs[2, i] + Y_t[i]
# }, 1:N)
# # Simpler: check if 0 is in the residual HPDI
# zero_in_hpdi <- sapply(1:N, function(i) {
#   0 >= resid_hpdi_per_obs[1, i] & 0 <= resid_hpdi_per_obs[2, i]
# })
# cat(sprintf("  89%% interval contains true value: %.1f%%\n", 100 * mean(zero_in_hpdi)))

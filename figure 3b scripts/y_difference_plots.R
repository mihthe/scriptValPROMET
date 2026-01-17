### accuracy of y

# run simulation file

# run analysis file

# library
library(rethinking)

# adapting proof of concept
Y_t <- sim$m_sur  # true mitotic counts on surgery
N <- length(Y_t)  # number of observations
l_s <- nrow(lambda)  # number of posterior samples 
Y_p <- matrix(nrow = N, ncol = l_s)  # predicted counts

# matrix with Poisson draws
for(i in 1:N) {
  Y_p[i, ] <- rpois(n = l_s, lambda = lambda[, i])
}

# compute residuals
resid <- matrix(nrow = N, ncol = l_s)

# residuals in matrix: true - predicted
for(i in 1:N) {
  resid[i, ] <- Y_t[i] - Y_p[i, ]
}

# HPDI
resid_hpdi <- HPDI(c(resid), prob = 0.89)

# Per-observation summaries
resid_median <- apply(resid, 1, median)
resid_hpdi_per_obs <- apply(resid, 1, HPDI, prob = 0.89)

# plot residual distribution
jpeg(paste0("output/figures/", "Residuals_distribution", ".jpg"), 
     units = "in", 
     width = 8, height = 6, res = 300)

hist(resid, freq = FALSE, 
     main = "Posterior Distribution of Residuals", 
     xlab = "Residual (True - Predicted)",
     col = scales::alpha("steelblue", 0.5),
     border = "white",
     breaks = 50)
abline(v = 0, col = "red", lwd = 2, lty = 2)
abline(v = resid_hpdi, col = "darkblue", lwd = 2, lty = 3)
legend("topright", 
       legend = c("Zero error", "89% HPDI"),
       col = c("red", "darkblue"),
       lwd = 2, lty = c(2, 3))

dev.off()

# median predicted values for each observation
lambda_median <- apply(lambda, 2, median)

# plot residuals vs fitted
jpeg(paste0("output/figures/", "Residuals_vs_Fitted", ".jpg"), 
     units = "in", 
     width = 10, height = 8, res = 300)
par(mfrow = c(2, 2))

plot(lambda_median, resid_median,
     xlab = "Fitted values (median lambda)",
     ylab = "Median residual",
     main = "Residuals vs Fitted Values",
     pch = 16, col = scales::alpha(1, 0.5))
abline(h = 0, col = "red", lwd = 2, lty = 2)
smooth_fit <- loess(resid_median ~ lambda_median)
lines(sort(lambda_median), predict(smooth_fit)[order(lambda_median)], 
      col = "blue", lwd = 2)

# plot absolute residuals vs fitted for heteroscedasticity
abs_resid_median <- apply(abs(resid), 1, median)
plot(lambda_median, abs_resid_median,
     xlab = "Fitted values (median lambda)",
     ylab = "Median absolute residual",
     main = "Scale-Location Plot",
     pch = 16, col = scales::alpha(1, 0.5))
smooth_fit2 <- loess(abs_resid_median ~ lambda_median)
lines(sort(lambda_median), predict(smooth_fit2)[order(lambda_median)], 
      col = "blue", lwd = 2)

# Q-Q plot of residuals
qqnorm(resid_median, 
       main = "Q-Q Plot of Median Residuals",
       pch = 16, col = scales::alpha(1, 0.5))
qqline(resid_median, col = "red", lwd = 2)

# histogram with normal overlay
hist(resid_median, freq = FALSE,
     main = "Distribution of Median Residuals",
     xlab = "Median Residual",
     col = scales::alpha("steelblue", 0.5),
     border = "white",
     breaks = 30)
curve(dnorm(x, mean(resid_median), sd(resid_median)),
      add = TRUE, col = "red", lwd = 2)
abline(v = 0, col = "darkred", lwd = 2, lty = 2)

dev.off()


# plot observed vs predicted
jpeg(paste0("output/figures/", "Observed_vs_Predicted", ".jpg"), 
     units = "in", 
     width = 10, height = 10, res = 300)
par(mfrow = c(1, 1))

# prediction intervals
pred_median <- apply(Y_p, 1, median)
pred_hpdi <- apply(Y_p, 1, HPDI, prob = 0.89)

plot(Y_t, pred_median,
     xlim = c(0, max(c(Y_t, pred_median)) + 2),
     ylim = c(0, max(c(Y_t, pred_median)) + 2),
     xlab = "Observed mitotic count (surgery)",
     ylab = "Predicted mitotic count (median)",
     main = "Observed vs Predicted Values\nwith 89% Prediction Intervals",
     pch = 16, col = scales::alpha(sim$L, 0.6),
     cex = 0.8)

# prediction intervals
segments(x0 = Y_t, 
         y0 = pred_hpdi[1, ], 
         y1 = pred_hpdi[2, ],
         col = scales::alpha(1, 0.2))

# perfect prediction line
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

# legend for locations
location_names <- c("Esophageal", "Gastric", "Duodenal", "Small intestine", "Colonrectum")
legend("topleft", 
       legend = location_names,
       pch = 16,
       col = scales::alpha(1:5, 0.6),
       title = "Location")

dev.off()


# plot residuals by covariates
jpeg(paste0("output/figures/", "Residuals_by_Covariates", ".jpg"), 
     units = "in", 
     width = 12, height = 8, res = 300)
par(mfrow = c(2, 3))

# plot residuals vs tumor size
plot(sim$Si, resid_median,
     xlab = "Tumor Size (standardized)",
     ylab = "Median Residual",
     main = "Residuals vs Tumor Size",
     pch = 16, col = scales::alpha(1, 0.5))
abline(h = 0, col = "red", lwd = 2, lty = 2)

# plot residuals vs biopsy surface
plot(sim$Su, resid_median,
     xlab = "Biopsy Surface (HPF)",
     ylab = "Median Residual",
     main = "Residuals vs Biopsy Surface",
     pch = 16, col = scales::alpha(1, 0.5))
abline(h = 0, col = "red", lwd = 2, lty = 2)

# plot residuals vs mitosis on biopsy
plot(sim$m_bio, resid_median,
     xlab = "Mitotic Count on Biopsy",
     ylab = "Median Residual",
     main = "Residuals vs Biopsy Count",
     pch = 16, col = scales::alpha(1, 0.5))
abline(h = 0, col = "red", lwd = 2, lty = 2)

# plot residuals by location 
boxplot(resid_median ~ sim$L,
        names = c("Esoph", "Stomach", "Duod", "Small", "Colon"),
        xlab = "Location",
        ylab = "Median Residual",
        main = "Residuals by Location",
        col = scales::alpha(2:5, 0.5))
abline(h = 0, col = "red", lwd = 2, lty = 2)

# plot residuals by NAC Response
boxplot(resid_median ~ sim$resp_Rx,
        names = c("No Response", "Response"),
        xlab = "NAC Response",
        ylab = "Median Residual",
        main = "Residuals by Treatment Response",
        col = scales::alpha(c(3, 4), 0.5))
abline(h = 0, col = "red", lwd = 2, lty = 2)

# plot residuals vs unknown confounder U
plot(sim$U, resid_median,
     xlab = "Unknown Confounder (U)",
     ylab = "Median Residual",
     main = "Residuals vs True Confounder",
     pch = 16, col = scales::alpha(1, 0.5))
abline(h = 0, col = "red", lwd = 2, lty = 2)
text(x = min(sim$U), y = max(resid_median),
     labels = "Note: U is unobserved in real data",
     pos = 4, col = "darkred", cex = 0.8)

dev.off()
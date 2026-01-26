## validation script for excel dataset

# run simlation file
# run analysis file
# save posteriors and standardisation of parametres
#saveRDS(post, "validation_output1/model_posterior.rds")
#saveRDS(std_params, "validation_output1/standardization_params.rds")
# run validation dataset file

# libraries
library(rethinking)
library(readxl)
library(ggplot2)
library(dplyr)
library(caret)
library(scales)  

# output directories if they don't exist
dir.create("validation_output1", showWarnings = FALSE)
dir.create("validation_output1/figures", showWarnings = FALSE)
dir.create("validation_output1/results", showWarnings = FALSE)

# =============================================================================
# PART 1: functions
# =============================================================================

# safe HPDI that handles constant values
safe_HPDI <- function(x, prob = 0.89) {
  if (length(unique(x)) == 1) {
    return(c(x[1], x[1]))
  }
  return(HPDI(x, prob))
}

# Miettinen & Lasota Risk Classification
classify_risk_ML <- function(size_mm, mitotic_count, location) {
  
  if(is.na(size_mm) | is.na(mitotic_count) | is.na(location)) {
    return(NA)
  }
  
  size_cm <- size_mm / 10
  is_gastric <- (location == 4)
  
  if(is_gastric) {
    if(mitotic_count <= 5) {
      if(size_cm <= 2) return("none")
      else if(size_cm <= 5) return("very_low")
      else if(size_cm <= 10) return("low")
      else return("moderate")
    } else {
      if(size_cm <= 2) return("none")
      else if(size_cm <= 5) return("moderate")
      else return("high")
    }
  } else {
    if(mitotic_count <= 5) {
      if(size_cm <= 2) return("none")
      else if(size_cm <= 5) return("low")
      else if(size_cm <= 10) return("moderate")
      else return("high")
    } else {
      return("high")
    }
  }
}

# =============================================================================
# PART 2: data loading and prediction functions
# =============================================================================

load_validation_data <- function(filepath = "validation_dataset_ready.xlsx") {
  
  validation_raw <- read_excel(filepath)
  validation_raw <- as.data.frame(validation_raw)
  
  required_cols <- c("patient_id", "size_mm", "L", "m_bio", 
                     "surface_mm2", "m_surg", "biopsy_surface_hpf")
  
  missing_cols <- setdiff(required_cols, names(validation_raw))
  if(length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s\nAvailable: %s", 
                 paste(missing_cols, collapse = ", "),
                 paste(names(validation_raw), collapse = ", ")))
  }
  
  std_params <- readRDS("validation_output1/standardization_params.rds")
  
  cat("  Standardization parameters loaded:\n")
  cat(sprintf("    Si_mean = %.3f, Si_sd = %.3f\n", std_params$Si_mean, std_params$Si_sd))
  cat(sprintf("    m_bio_mean = %.3f, m_bio_sd = %.3f\n", std_params$m_bio_mean, std_params$m_bio_sd))
  cat(sprintf("    Su_mean = %.3f, Su_sd = %.3f\n", std_params$Su_mean, std_params$Su_sd))
  
  validation_data <- validation_raw %>%
    mutate(
      size_std = (size_mm - std_params$Si_mean) / std_params$Si_sd,
      biopsy_mitosis_raw = m_bio,
      biopsy_surface_raw = biopsy_surface_hpf,
      data_source = "external_validation"
    )

  return(list(
    data = validation_data,
    std_params = std_params
  ))
}

predict_mitosis <- function(patient_data, posterior, std_params) {
  
  m_bio_std <- (patient_data$m_bio - std_params$m_bio_mean) / std_params$m_bio_sd
  Su_std <- (patient_data$biopsy_surface_hpf - std_params$Su_mean) / std_params$Su_sd
  Si_std <- patient_data$size_std
  L <- patient_data$L
  
  n_samples <- length(posterior$a)
  log_lambda <- numeric(n_samples)
  
  for(i in 1:n_samples) {
    log_lambda[i] <- posterior$a[i] + 
      posterior$b[i, L] * Si_std +
      posterior$g[i] * m_bio_std +
      posterior$e[i, L] * Su_std
  }
  
  lambda <- exp(log_lambda)
  predicted_counts <- rpois(n_samples, lambda = lambda)
  
  return(predicted_counts)
}

predict_all_validation <- function(validation_data, posterior, std_params) {
  
  n_patients <- nrow(validation_data)
  n_samples <- length(posterior$a)
  
  cat("MAKING PREDICTIONS\n")
  cat(sprintf("  Patients: %d\n", n_patients))
  cat(sprintf("  Posterior samples: %d\n", n_samples))
  
  predicted_samples <- matrix(NA, nrow = n_samples, ncol = n_patients)
  
  for(i in 1:n_patients) {
    predicted_samples[, i] <- predict_mitosis(
      patient_data = validation_data[i, ],
      posterior = posterior,
      std_params = std_params
    )
    
    if(i %% 10 == 0 | i == n_patients) {
      cat(sprintf("\r  Progress: %d/%d", i, n_patients))
    }
  }
  
  predicted_mean <- colMeans(predicted_samples)
  predicted_median <- apply(predicted_samples, 2, median)
  
  predicted_ci_lower <- numeric(n_patients)
  predicted_ci_upper <- numeric(n_patients)
  
  for(i in 1:n_patients) {
    hpdi_result <- safe_HPDI(predicted_samples[, i], 0.89)
    predicted_ci_lower[i] <- hpdi_result[1]
    predicted_ci_upper[i] <- hpdi_result[2]
  }
  
  validation_data <- validation_data %>%
    mutate(
      predicted_mean = predicted_mean,
      predicted_median = predicted_median,
      predicted_ci_lower = predicted_ci_lower,
      predicted_ci_upper = predicted_ci_upper,
      prediction_error = predicted_mean - m_surg,
      abs_error = abs(prediction_error),
      in_ci = (m_surg >= predicted_ci_lower) & (m_surg <= predicted_ci_upper)
    )
  
  return(list(
    data = validation_data,
    posterior_samples = predicted_samples
  ))
}

# =============================================================================
# PART 3: metrics
# =============================================================================

calculate_validation_metrics <- function(validation_results) {
  
  data <- validation_results$data
  
  cat("CALCULATING VALIDATION METRICS\n")
  cat(paste(rep("─", 50), collapse = ""), "\n")
  
  mae <- mean(abs(data$prediction_error))
  rmse <- sqrt(mean(data$prediction_error^2))
  bias <- mean(data$prediction_error)
  r_squared <- cor(data$predicted_mean, data$m_surg)^2
  coverage_89 <- mean(data$in_ci)
  
  risk_levels <- c("none", "very_low", "low", "moderate", "high")
  
  actual_risk <- mapply(classify_risk_ML,
                        size_mm = data$size_mm,
                        mitotic_count = data$m_surg,
                        location = data$L)
  
  predicted_risk <- mapply(classify_risk_ML,
                           size_mm = data$size_mm,
                           mitotic_count = round(data$predicted_mean),
                           location = data$L)
  
  biopsy_risk <- mapply(classify_risk_ML,
                        size_mm = data$size_mm,
                        mitotic_count = data$m_bio,
                        location = data$L)
  
  actual_risk <- factor(actual_risk, levels = risk_levels)
  predicted_risk <- factor(predicted_risk, levels = risk_levels)
  biopsy_risk <- factor(biopsy_risk, levels = risk_levels)
  
  cm_model <- confusionMatrix(predicted_risk, actual_risk)
  cm_biopsy <- confusionMatrix(biopsy_risk, actual_risk)
  
  metrics <- list(
    continuous = list(
      mae = mae,
      rmse = rmse,
      bias = bias,
      r_squared = r_squared,
      coverage_89 = coverage_89,
      n = nrow(data)
    ),
    
    classification = list(
      model_accuracy = cm_model$overall["Accuracy"],
      biopsy_accuracy = cm_biopsy$overall["Accuracy"],
      model_kappa = cm_model$overall["Kappa"],
      biopsy_kappa = cm_biopsy$overall["Kappa"],
      cm_model = cm_model$table,
      cm_biopsy = cm_biopsy$table
    ),
    
    risks = list(
      actual = actual_risk,
      predicted = predicted_risk,
      biopsy = biopsy_risk
    )
  )
  
  cat("\n")
  cat("CONTINUOUS METRICS:\n")
  cat(sprintf("  Mean Absolute Error (MAE):     %.2f mitoses\n", mae))
  cat(sprintf("  Root Mean Square Error (RMSE): %.2f mitoses\n", rmse))
  cat(sprintf("  Bias (mean error):             %.2f mitoses\n", bias))
  cat(sprintf("  R² (correlation squared):      %.3f\n", r_squared))
  cat(sprintf("  89%% HPDI Coverage:             %.1f%%\n", coverage_89 * 100))
  
  cat("\nCLASSIFICATION METRICS (Miettinen & Lasota):\n")
  cat(sprintf("  PROMETheus Model Accuracy: %.1f%% (κ = %.3f)\n", 
              cm_model$overall["Accuracy"] * 100, cm_model$overall["Kappa"]))
  cat(sprintf("  Biopsy Alone Accuracy:     %.1f%% (κ = %.3f)\n", 
              cm_biopsy$overall["Accuracy"] * 100, cm_biopsy$overall["Kappa"]))
  
  improvement <- (cm_model$overall["Accuracy"] - cm_biopsy$overall["Accuracy"]) * 100
  cat(sprintf("\n  ➜ Model improvement over biopsy: %+.1f percentage points\n", improvement))
  
  cat("\n")
  
  return(metrics)
}

# =============================================================================
# PART 4: plots
# =============================================================================

# location colors 
loc_colors <- c("coral", "gold", "forestgreen", "steelblue")
loc_names <- c("Colon-rectum", "Duodenum", "Small intestine", "Stomach")

# -----------------------------------------------------------------------------
# Plot 1: Calibration Plot
# -----------------------------------------------------------------------------
create_calibration_plot <- function(validation_results, metrics, 
                                    output_dir = "validation_output1/figures") {
  
  data <- validation_results$data
  point_colors <- loc_colors[data$L]
  max_val <- max(c(data$predicted_mean, data$m_surg), na.rm = TRUE)
  
  jpeg(file.path(output_dir, "calibration_plot.jpg"),
       width = 10, height = 9, units = "in", res = 300)
  
  par(mar = c(8, 5, 4, 2))
  
  plot(data$predicted_mean, data$m_surg,
       xlab = "Predicted Mitotic Count (PROMETheus)",
       ylab = "Observed Mitotic Count (Surgery)",
       main = "External Validation: Calibration Plot",
       pch = 16, col = scales::alpha(point_colors, 0.7), cex = 1.8,
       xlim = c(0, max_val * 1.1), ylim = c(0, max_val * 1.1))
  
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  
  segments(x0 = data$predicted_ci_lower,
           y0 = data$m_surg,
           x1 = data$predicted_ci_upper,
           y1 = data$m_surg,
           col = rgb(0, 0, 0, 0.3), lwd = 1.5)
  
  # legend("topleft",
  #        legend = c(
  #          sprintf("N = %d", metrics$continuous$n),
  #          sprintf("MAE = %.2f", metrics$continuous$mae),
  #          sprintf("RMSE = %.2f", metrics$continuous$rmse),
  #          sprintf("R² = %.3f", metrics$continuous$r_squared),
  #          sprintf("89%% CI Coverage = %.1f%%", metrics$continuous$coverage_89 * 100)
  #        ),
  #        bty = "n", cex = 1.0)
  
  par(xpd = TRUE)
  legend("bottom",
         inset = c(0, -0.25),
         legend = c("Perfect prediction", loc_names),
         lty = c(2, rep(0, 4)),
         lwd = c(2, rep(0, 4)),
         pch = c(NA, rep(16, 4)),
         col = c("red", loc_colors),
         pt.cex = c(1, rep(1.8, 4)),
         bty = "n",
         horiz = TRUE,
         cex = 0.95)
  par(xpd = FALSE)
  
  dev.off()
}

# -----------------------------------------------------------------------------
# Plot 2: Residuals Plot
# -----------------------------------------------------------------------------
create_residuals_plot <- function(validation_results, 
                                  output_dir = "validation_output1/figures") {
  
  data <- validation_results$data
  point_colors <- loc_colors[data$L]
  
  jpeg(file.path(output_dir, "residuals_plot.jpg"),
       width = 10, height = 8, units = "in", res = 300)
  
  par(mar = c(7, 5, 4, 2))
  
  plot(data$predicted_mean, data$prediction_error,
       xlab = "Predicted Mitotic Count",
       ylab = "Residual (Predicted - Observed)",
       main = "Residual Plot",
       pch = 16, col = scales::alpha(point_colors, 0.7), cex = 1.5)
  
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  abline(h = c(-5, 5), col = "orange", lwd = 1, lty = 3)
  
  par(xpd = TRUE)
  legend("bottom",
         inset = c(0, -0.22),
         legend = c("Zero error", "±5 mitoses", loc_names),
         lty = c(2, 3, rep(0, 4)),
         lwd = c(2, 1, rep(0, 4)),
         pch = c(NA, NA, rep(16, 4)),
         col = c("red", "orange", loc_colors),
         pt.cex = c(1, 1, rep(1.5, 4)),
         bty = "n",
         horiz = TRUE,
         cex = 0.9)
  par(xpd = FALSE)
  
  dev.off()
  cat("✓ Residuals plot saved (JPG)\n")
}

# -----------------------------------------------------------------------------
# Plot 3: Posterior Predictive Check 
# -----------------------------------------------------------------------------
create_posterior_predictive_plot <- function(validation_results, 
                                             output_dir = "validation_output1/figures",
                                             y_cap = 50) {  # Cap y-axis at 50
  
  data <- validation_results$data
  n <- nrow(data)
  
  j <- order(data$predicted_median)
  loc_colors <- c("coral", "gold", "forestgreen", "steelblue")
  loc_names <- c("Colon-rectum", "Duodenum", "Small intestine", "Stomach")
  loc_colors_bg <- scales::alpha(loc_colors[data$L[j]], 0.15)
  
  # identify cases that exceed the cap
  exceeds_cap <- data$predicted_ci_upper[j] > y_cap | data$m_surg[j] > y_cap
  n_exceeds <- sum(exceeds_cap)
  
  jpeg(file.path(output_dir, "posterior_predictive_check.jpg"),
       width = 14, height = 7, units = "in", res = 300)
  
  par(mar = c(7, 5, 4, 2))
  
  # use capped y-axis
  plot(NULL, 
       xlim = c(0.5, n + 0.5),
       ylim = c(-0.5, y_cap),
       xlab = "",
       ylab = "Mitotic Count (per 5 mm²)",
       main = "Posterior Predictive Check - External Validation",
       xaxt = "n",
       bty = "n")
   
  axis(1, at = 1:n, labels = 1:n, las = 2, cex.axis = 0.7)
  mtext("Cases (ordered by predicted count)", side = 1, line = 3.5)
  
  # background strips for location
  for(i in 1:n) {
    rect(i - 0.4, -1, i + 0.4, y_cap + 5, 
         col = loc_colors_bg[i], border = NA)
  }
  
  # risk threshold line
  abline(h = 5, col = "gray50", lty = 2, lwd = 1)
  text(0, 50, "Risk threshold (5 mitoses)", adj = 0, cex = 0.8, col = "gray40", font = 3)
  
  # cap the CI values for plotting
  ci_lower_plot <- pmin(data$predicted_ci_lower[j], y_cap)
  ci_upper_plot <- pmin(data$predicted_ci_upper[j], y_cap)
  median_plot <- pmin(data$predicted_median[j], y_cap)
  surg_plot <- pmin(data$m_surg[j], y_cap)
  bio_plot <- pmin(data$m_bio[j], y_cap)
  
  # 89% CI segments
  segments(x0 = 1:n, y0 = ci_lower_plot,
           y1 = ci_upper_plot, lwd = 2.5, col = "gray40")
  
  # arrows for CIs that exceed the cap
  for(i in 1:n) {
    if(data$predicted_ci_upper[j][i] > y_cap) {
      arrows(i, y_cap - 2, i, y_cap, length = 0.08, col = "gray40", lwd = 2)
    }
  }
  
  # predicted median
  points(1:n, median_plot, pch = 16, cex = 1.3, col = "black")
  
  # observed surgery count
  points(1:n, surg_plot, pch = 4, cex = 1.3, col = "red", lwd = 2)
  
  # biopsy count
  points(1:n, bio_plot, pch = 1, cex = 1, col = "blue", lwd = 1.5)
  
  # note cases that exceed the cap
  if(n_exceeds > 0) {
    exceed_indices <- which(exceeds_cap)
    for(idx in exceed_indices) {
      # text annotation showing actual values
      actual_median <- round(data$predicted_median[j][idx])
      actual_surg <- data$m_surg[j][idx]
      text(idx, y_cap - 3, 
           sprintf("λ=%d\ny=%d", actual_median, actual_surg),
           cex = 0.6, col = "darkred", font = 2)
    }
  }
  
  # legend below plot
  par(xpd = TRUE)
  legend("bottom", 
         inset = c(0, -0.28),
         legend = c("Predicted (median)", "89% HPDI", "Observed (surgery)", "Biopsy count",
                    loc_names),
         pch = c(16, NA, 4, 1, 22, 22, 22, 22),
         lty = c(0, 1, 0, 0, 0, 0, 0, 0),
         lwd = c(0, 2.5, 2, 1.5, 0, 0, 0, 0),
         col = c("black", "gray40", "red", "blue", "black", "black", "black", "black"),
         pt.bg = c(NA, NA, NA, NA, scales::alpha(loc_colors, 0.4)),
         pt.cex = c(1.3, 1, 1.3, 1, 2, 2, 2, 2),
         bty = "n", 
         horiz = TRUE,
         cex = 0.85)
  par(xpd = FALSE)
  
  dev.off()
  
  # report
  cat("✓ Posterior predictive check plot saved (JPG)\n")
  if(n_exceeds > 0) {
    cat(sprintf("  Note: %d case(s) exceed y-axis cap of %d\n", n_exceeds, y_cap))
  }
}

# # -----------------------------------------------------------------------------
# # Plot 3: Posterior Predictive Check - larger y-axis and outlier case 44
# # -----------------------------------------------------------------------------
# create_posterior_predictive_plot <- function(validation_results, 
#                                              output_dir = "validation_output1/figures") {
#   
#   data <- validation_results$data
#   n <- nrow(data)
#   
#   j <- order(data$predicted_median)
#   loc_colors_bg <- scales::alpha(loc_colors[data$L[j]], 0.15)
#   y_max <- max(c(data$m_surg, data$predicted_ci_upper), na.rm = TRUE) + 3
#   
#   jpeg(file.path(output_dir, "posterior_predictive_check.jpg"),
#        width = 14, height = 7, units = "in", res = 300)
#   
#   par(mar = c(7, 5, 4, 2))
#   
#   # x-axis properly aligned to n (number of cases)
#   plot(NULL, 
#        xlim = c(0.5, n + 0.5),
#        ylim = c(-0.5, y_max),
#        xlab = "",
#        ylab = "Mitotic Count (per 5 mm²)",
#        main = "Posterior Predictive Check - External Validation",
#        xaxt = "n",
#        bty = "n")
#   
#   axis(1, at = 1:n, labels = 1:n, las = 2, cex.axis = 0.7)
#   mtext("Cases (ordered by predicted count)", side = 1, line = 3.5)
#   
#   for(i in 1:n) {
#     rect(i - 0.4, -1, i + 0.4, y_max + 5, 
#          col = loc_colors_bg[i], border = NA)
#   }
#   
#   # rRisk threshold text sx
#   abline(h = 5, col = "gray50", lty = 2, lwd = 1.5)
#   mtext("- - - Risk threshold (5 mitoses)", side = 3, line = -1.5, 
#         adj = 0, cex = 0.8, col = "gray40", font = 3)
#   
#   segments(x0 = 1:n, y0 = data$predicted_ci_lower[j],
#            y1 = data$predicted_ci_upper[j], lwd = 2.5, col = "gray40")
#   
#   points(1:n, data$predicted_median[j], pch = 16, cex = 1.3, col = "black")
#   points(1:n, data$m_surg[j], pch = 4, cex = 1.3, col = "red", lwd = 2)
#   points(1:n, data$m_bio[j], pch = 1, cex = 1, col = "blue", lwd = 1.5)
#   
#   # legend below plot
#   par(xpd = TRUE)
#   legend("bottom", 
#          inset = c(0, -0.28),
#          legend = c("Predicted (median)", "89% HPDI", "Observed (surgery)", "Biopsy count",
#          loc_names),
#          pch = c(16, NA, 4, 1, 22, 22, 22, 22),
#          lty = c(0, 1, 0, 0, 0, 0, 0, 0),
#          lwd = c(0, 2.5, 2, 1.5, 0, 0, 0, 0),
#          col = c("black", "gray40", "red", "blue", "black", "black", "black", "black"),
#          pt.bg = c(NA, NA, NA, NA, scales::alpha(loc_colors, 0.4)),
#          pt.cex = c(1.3, 1, 1.3, 1, 2, 2, 2, 2),
#          bty = "n", 
#          horiz = TRUE,
#          cex = 0.85)
#   par(xpd = FALSE)
#   
#   dev.off()
# }

# # -----------------------------------------------------------------------------
# # Plot 4: Prediction Error by Tumor Size (CORRECTED)
# # -----------------------------------------------------------------------------
# create_error_by_size_plot <- function(validation_results, 
#                                       output_dir = "validation_output1/figures") {
#   
#   data <- validation_results$data
#   n <- nrow(data)
#   
#   pred_error <- data$predicted_median - data$m_surg
#   j <- order(data$size_mm)
#   
#   loc_colors_bg <- scales::alpha(loc_colors[data$L[j]], 0.2)
#   is_over <- pred_error[j] > 0
#   line_colors <- ifelse(is_over, "forestgreen", "red")
#   
#   mae <- mean(abs(pred_error))
#   bias <- mean(pred_error)
#   n_over <- sum(pred_error > 0)
#   n_under <- sum(pred_error < 0)
#   
#   y_lim <- max(abs(pred_error), na.rm = TRUE) * 1.1
#   
#   jpeg(file.path(output_dir, "prediction_error_by_size.jpg"),
#        width = 12, height = 8, units = "in", res = 300)
#   
#   # margin  
#   par(mar = c(10, 5, 4, 2))
#   
#   plot(NULL,
#        xlim = c(0.5, n + 0.5),
#        ylim = c(-y_lim, y_lim * 0.15),
#        xlab = "",
#        ylab = expression("Prediction Error (Predicted " * lambda * " - Observed Mitotic Count)"),
#        main = "External Validation: Prediction Error Ordered by Tumor Size",
#        xaxt = "n",
#        bty = "n")
#   
#   for(i in 1:n) {
#     rect(i - 0.4, -y_lim * 1.5, i + 0.4, y_lim * 0.5,
#          col = loc_colors_bg[i], border = NA)
#   }
#   
#   abline(h = 0, col = "black", lwd = 2)
#   
#   segments(x0 = 1:n, y0 = 0, y1 = pred_error[j],
#            col = line_colors, lwd = 3)
#   
#   axis(1, at = 1:n, labels = round(data$size_mm[j]), las = 2, cex.axis = 0.7)
#   mtext("Tumor Size (mm)", side = 1, line = 4)
#   
#   # summary statistics box
#   text_x <- 1
#   text_y <- y_lim * 0.1
#   text(text_x, text_y, pos = 4, cex = 0.9, font = 2, labels = "Summary Statistics")
#   text(text_x, text_y - y_lim * 0.06, pos = 4, cex = 0.85,
#        labels = sprintf("N = %d patients", n))
#   text(text_x, text_y - y_lim * 0.11, pos = 4, cex = 0.85,
#        labels = sprintf("MAE = %.2f mitoses", mae))
#   text(text_x, text_y - y_lim * 0.16, pos = 4, cex = 0.85,
#        labels = sprintf("Bias = %.2f mitoses", bias))
#   text(text_x, text_y - y_lim * 0.21, pos = 4, cex = 0.85,
#        labels = sprintf("Overestimated: %d (%.1f%%)", n_over, 100*n_over/n))
#   text(text_x, text_y - y_lim * 0.26, pos = 4, cex = 0.85,
#        labels = sprintf("Underestimated: %d (%.1f%%)", n_under, 100*n_under/n))
#   
#   # legends below plot
#   par(xpd = TRUE)
#   
#   legend("bottom",
#          inset = c(0, -0.35),
#          title = "Tumor Location",
#          legend = loc_names,
#          pch = 22,
#          pt.bg = scales::alpha(loc_colors, 0.4),
#          col = "black",
#          pt.cex = 2.5,
#          bty = "n",
#          horiz = TRUE,
#          cex = 0.9)
#   
#   legend("bottom",
#          inset = c(0, -0.45),
#          title = "Prediction Direction",
#          legend = c(expression("Overestimation (" * lambda > y * ")"),
#                     expression("Underestimation (" * lambda < y * ")"),
#                     "Perfect prediction"),
#          lty = c(1, 1, 1),
#          lwd = c(3, 3, 2),
#          col = c("forestgreen", "red", "black"),
#          bty = "n",
#          horiz = TRUE,
#          cex = 0.9)
#   
#   par(xpd = FALSE)
#   
#   dev.off()
# }

# -----------------------------------------------------------------------------
# Plot 4: Prediction Error by Tumor Size 
# -----------------------------------------------------------------------------
create_error_by_size_plot <- function(validation_results, 
                                      output_dir = "validation_output1/figures") {
  
  data <- validation_results$data
  n <- nrow(data)
  
  pred_error <- data$predicted_median - data$m_surg
  j <- order(data$size_mm)
  
  loc_colors <- c("coral", "gold", "forestgreen", "steelblue")
  loc_names <- c("Colon-rectum", "Duodenum", "Small intestine", "Stomach")
  loc_colors_bg <- scales::alpha(loc_colors[data$L[j]], 0.2)
  
  is_over <- pred_error[j] > 0
  line_colors <- ifelse(is_over, "forestgreen", "red")
  
  mae <- mean(abs(pred_error))
  bias <- mean(pred_error)
  n_over <- sum(pred_error > 0)
  n_under <- sum(pred_error < 0)
  n_exact <- sum(pred_error == 0)
  
  # y-axis SYMMETRIC around zero
  y_max_abs <- max(abs(pred_error), na.rm = TRUE) * 1.15
  
  # figure 
  jpeg(file.path(output_dir, "prediction_error_by_size.jpg"),
       width = 12, height = 10, units = "in", res = 300)  
  
  # margin 
  par(mar = c(12, 5, 4, 2))  
  
  # plot
  plot(NULL,
       xlim = c(0.5, n + 0.5),
       ylim = c(-y_max_abs, y_max_abs),  
       xlab = "",
       ylab = expression("Prediction Error (Predicted " * lambda * " - Observed Mitotic Count)"),
       main = "External Validation: Prediction Error Ordered by Tumor Size",
       xaxt = "n",
       bty = "n")
  
  # background strips for location
  for(i in 1:n) {
    rect(i - 0.4, -y_max_abs * 1.5, i + 0.4, y_max_abs * 1.5,
         col = loc_colors_bg[i], border = NA)
  }
  
  # zero line (perfect prediction)
  abline(h = 0, col = "black", lwd = 2)
  
  # reference lines for context
  abline(h = c(-10, 10), col = "gray70", lty = 3, lwd = 1)
  abline(h = c(-5, 5), col = "gray50", lty = 2, lwd = 1)
  
  # error segments
  segments(x0 = 1:n, y0 = 0, y1 = pred_error[j],
           col = line_colors, lwd = 3)
  
  # x-axis: tumor sizes
  axis(1, at = 1:n, labels = round(data$size_mm[j]), las = 2, cex.axis = 0.7)
  mtext("Tumor Size (mm)", side = 1, line = 4)
  
  # # Summary statistics box (positioned in upper right, away from data)
  # text_x <- n * 0.75
  # text_y <- y_max_abs * 0.85
  # 
  # # Add a semi-transparent background box
  # rect(text_x - 2, text_y - y_max_abs * 0.35, n + 0.5, text_y + y_max_abs * 0.1,
  #      col = rgb(1, 1, 1, 0.8), border = NA)
  # 
  # text(text_x, text_y, pos = 4, cex = 0.9, font = 2, labels = "Summary Statistics")
  # text(text_x, text_y - y_max_abs * 0.07, pos = 4, cex = 0.85,
  #      labels = sprintf("N = %d patients", n))
  # text(text_x, text_y - y_max_abs * 0.13, pos = 4, cex = 0.85,
  #      labels = sprintf("MAE = %.2f mitoses", mae))
  # text(text_x, text_y - y_max_abs * 0.19, pos = 4, cex = 0.85,
  #      labels = sprintf("Bias = %.2f mitoses", bias))
  # text(text_x, text_y - y_max_abs * 0.25, pos = 4, cex = 0.85,
  #      labels = sprintf("Overestimated: %d (%.1f%%)", n_over, 100*n_over/n))
  # text(text_x, text_y - y_max_abs * 0.31, pos = 4, cex = 0.85,
  #      labels = sprintf("Underestimated: %d (%.1f%%)", n_under, 100*n_under/n))
  
  # legends positioned with proper spacing
  par(xpd = TRUE)
  
  # first legend row: Tumor Location
  legend("bottom",
         inset = c(0, -0.22),  
         title = "Tumor Location",
         title.font = 2,
         legend = loc_names,
         pch = 22,
         pt.bg = scales::alpha(loc_colors, 0.4),
         col = "black",
         pt.cex = 2.5,
         bty = "n",
         horiz = TRUE,
         cex = 0.9)
  
  # second legend row: Prediction Direction 
  legend("bottom",
         inset = c(0, -0.32),  
         title = "Prediction Direction",
         title.font = 2,
         legend = c("Overestimation (λ > y)",
                    "Underestimation (λ < y)",
                    "Perfect prediction"),
         lty = c(1, 1, 1),
         lwd = c(3, 3, 2),
         col = c("forestgreen", "red", "black"),
         bty = "n",
         horiz = TRUE,
         cex = 0.9)
  
  par(xpd = FALSE)
  
  dev.off()
}

# -----------------------------------------------------------------------------
# Plot 5: Risk Classification Heatmaps
# -----------------------------------------------------------------------------
create_risk_classification_plot <- function(metrics, 
                                            output_dir = "validation_output1/figures") {
  
  jpeg(file.path(output_dir, "risk_classification.jpg"),
       width = 14, height = 7, units = "in", res = 300)
  
  par(mfrow = c(1, 2), mar = c(8, 8, 4, 2))
  
  plot_confusion_heatmap <- function(cm, title, color_palette) {
    
    row_sums <- rowSums(cm)
    col_sums <- colSums(cm)
    keep <- (row_sums > 0) | (col_sums > 0)
    cm_sub <- cm[keep, keep, drop = FALSE]
    
    if(nrow(cm_sub) == 0) {
      plot.new()
      text(0.5, 0.5, "No data", cex = 2)
      return()
    }
    
    cm_prop <- prop.table(cm_sub, margin = 2) * 100
    cm_prop[is.nan(cm_prop)] <- 0
    
    n_classes <- nrow(cm_sub)
    
    image(1:n_classes, 1:n_classes, 
          t(cm_prop)[, n_classes:1, drop = FALSE],
          col = colorRampPalette(color_palette)(50),
          xlab = "", ylab = "",
          main = title,
          axes = FALSE)
    
    axis(1, at = 1:n_classes, labels = rownames(cm_sub), las = 2, cex.axis = 1)
    axis(2, at = 1:n_classes, labels = rev(colnames(cm_sub)), las = 1, cex.axis = 1)
    mtext("Predicted Risk", side = 1, line = 5, cex = 1.1)
    mtext("Actual Risk (Surgery)", side = 2, line = 5, cex = 1.1)
    
    abline(h = seq(0.5, n_classes + 0.5, 1), col = "white", lwd = 0.5)
    abline(v = seq(0.5, n_classes + 0.5, 1), col = "white", lwd = 0.5)
    
    for(i in 1:n_classes) {
      for(j in 1:n_classes) {
        val <- cm_sub[i, n_classes - j + 1]
        pct <- cm_prop[i, n_classes - j + 1]
        text_col <- ifelse(pct > 40, "white", "black")
        text(i, j, sprintf("%d\n(%.0f%%)", val, pct),
             col = text_col, cex = 1.1, font = 2)
      }
    }
  }
  
  plot_confusion_heatmap(metrics$classification$cm_model, 
                         "PROMETheus Model", 
                         c("white", "steelblue", "darkblue"))
  
  plot_confusion_heatmap(metrics$classification$cm_biopsy, 
                         "Biopsy Alone (Current Practice)", 
                         c("white", "coral", "darkred"))
  
  dev.off()
  cat("✓ Risk classification plot saved (JPG)\n")
}

# -----------------------------------------------------------------------------
# create all plots function
# -----------------------------------------------------------------------------
create_all_validation_plots <- function(validation_results, metrics, 
                                        output_dir = "validation_output1/figures") {
  
  cat("\nCREATING VALIDATION PLOTS (JPG format)\n")
  cat(paste(rep("─", 50), collapse = ""), "\n")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  create_calibration_plot(validation_results, metrics, output_dir)
  create_residuals_plot(validation_results, output_dir)
  create_posterior_predictive_plot(validation_results, output_dir)
  create_error_by_size_plot(validation_results, output_dir)
  create_risk_classification_plot(metrics, output_dir)
  
  cat("\n✓ All plots saved to:", output_dir, "\n\n")
}

# =============================================================================
# PART 5: MAIN VALIDATION WORKFLOW
# =============================================================================

run_validation <- function(validation_file = "validation_dataset_ready.xlsx") {
  
  cat("\n")
  cat("PROMETheus EXTERNAL VALIDATION")

  # load trained model posterior
  posterior <- readRDS("validation_output1/model_posterior.rds")
  cat("✓ Loaded trained model posterior\n")
  cat(sprintf("  Posterior samples: %d\n", length(posterior$a)))
  
  # load validation data
  loaded <- load_validation_data(validation_file)
  validation_data <- loaded$data
  std_params <- loaded$std_params
  
  # make predictions
  validation_results <- predict_all_validation(validation_data, posterior, std_params)
  
  # calculate metrics
  metrics <- calculate_validation_metrics(validation_results)
  
  # create ALL plots 
  create_all_validation_plots(validation_results, metrics)
  
  # save results
  results <- list(
    data = validation_results$data,
    posterior_samples = validation_results$posterior_samples,
    metrics = metrics,
    info = list(
      n_patients = nrow(validation_data),
      validation_file = validation_file,
      date = Sys.Date()
    )
  )
  
  saveRDS(results, "validation_output1/results/validation_results.rds")
  write.csv(validation_results$data, 
            "validation_output1/results/predictions.csv", 
            row.names = FALSE)
  
  return(results)
}

# =============================================================================
# RUN THE VALIDATION
# =============================================================================

# runs the complete validation and creates all plots
results <- run_validation("validation_dataset_ready.xlsx")

# The 'results' object now contains everything:
#   results$data              - validation data with predictions
#   results$posterior_samples - all posterior samples
#   results$metrics           - all calculated metrics
#   results$info              - metadata

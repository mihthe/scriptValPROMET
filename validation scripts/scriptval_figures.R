# figures for Results section using real validation data
################################################################################

#libraries
library(scales)
library(ggplot2)
library(dplyr)

################################################################################
# FIRST run
# data simulation file to train model
# data analysis file to fit model
# validation file to load results
################################################################################

################################################################################
# load validation results
################################################################################

if(!exists("results")) {
  if(file.exists("validation_output1/results/validation_results.rds")) {
    results <- readRDS("validation_output1/results/validation_results.rds")
    cat("✓ Loaded validation results\n")
  } else {
    stop("Run validation first: results <- run_validation()")
  }
}

validation_data <- results$data 
predicted_samples <- results$posterior_samples 

################################################################################
# FIGURE 1: Posterior Predictive Check (like fig1.R but for validation)
################################################################################

create_figure1_validation <- function(validation_data, predicted_samples) {
  
  cat("\nCreating Figure 1: Posterior Predictive Check...\n")
  
  N <- nrow(validation_data)
  
  # calculate median predictions
  lambda_median <- apply(predicted_samples, 2, median)
  
  # order patients by median prediction
  idx <- order(lambda_median)
  
  # Y-axis limits
  y_min <- 0
  y_max <- max(c(validation_data$surgery_mitosis, 
                 apply(predicted_samples, 2, quantile, probs = 0.975)))
  
  # plot
  jpeg("validation_output1/figures/figure1_validation.jpg", 
       units = "in", width = 14, height = 8, res = 300)
  
  par(mai = c(0.9, 0.9, 0.9, 0.4))
  
  plot(NULL, 
       xlim = c(0.5, N + 0.5),  
       ylim = c(y_min, y_max),  
       xlab = 'Validation Patient (ordered by predicted mitotic count)', 
       ylab = 'Mitotic count (per 5 mm²)',
       main = 'Figure 1: Posterior Predictive Check - External Validation',
       bty = 'l',
       cex.lab = 1.2,
       cex.main = 1.4)
  
  # reference line at mitotic count = 5
  abline(h = 5, lty = 2, col = "gray40", lwd = 2)
  text(N * 0.95, 5.5, "High-risk threshold", pos = 3, cex = 0.9, col = "gray40")
  
  # scaling factor for violin width
  s_factor <- 0.4
  
  # for each patient
  for(i in 1:N) {
    patient_idx <- idx[i]
    
    # posterior distribution for this patient
    patient_predictions <- predicted_samples[, patient_idx]
    
    # density of posterior
    dens_obj <- density(patient_predictions, adjust = 1.2)
    x_dens <- dens_obj$y
    y_dens <- dens_obj$x
    
    # normalize density
    x_dens_norm <- x_dens / max(x_dens) * s_factor
    
    # dx half of violin
    polygon(i + x_dens_norm, y_dens, 
            col = alpha("skyblue", 0.5), 
            border = alpha("steelblue", 0.7), 
            lwd = 0.3)
    
    # sx half of violin
    polygon(i - x_dens_norm, y_dens, 
            col = alpha("skyblue", 0.5), 
            border = alpha("steelblue", 0.7), 
            lwd = 0.3)
    
    # observed mitotic count on surgery
    points(i, validation_data$surgery_mitosis[patient_idx], 
           pch = 4, col = "red", lwd = 2.5, cex = 1.3)
    
    # median prediction
    segments(i - 0.35, lambda_median[patient_idx], 
             i + 0.35, lambda_median[patient_idx], 
             lwd = 2.5, col = "black")
  }
  
  # legend
  legend('topleft', 
         pch = c(4, 15, NA),
         lty = c(NA, NA, 1),
         col = c("red", "lightblue", "black"),
         pt.cex = c(1.3, 2, NA),
         lwd = c(NA, NA, 2.5),
         legend = c("Observed (surgery specimen)", 
                    "Posterior predictive distribution",
                    "Median prediction"),
         bty = "n",
         cex = 1.1)
  
  # additional info
  mtext(sprintf("N = %d patients | External validation cohort", N), 
        side = 3, line = 0.5, cex = 0.9, adj = 1)
  
  dev.off()
  
  cat("✓ Figure 1 saved: validation_output1/figures/figure1_validation.jpg\n")
}

################################################################################
# FIGURE 2A: Prediction Error by Biopsy Surface
################################################################################

create_figure2a_validation <- function(validation_data, predicted_samples) {
  
  cat("\nCreating Figure 2A: Prediction Error by Biopsy Surface...\n")
  
  # calculate differences
  difference <- validation_data$predicted_median - validation_data$surgery_mitosis
  abs_difference <- abs(difference)
  
  # create comprehensive data frame
  plot_data <- data.frame(
    patient_id = 1:nrow(validation_data),
    difference = difference,
    abs_difference = abs_difference,
    surface = validation_data$biopsy_surface_hpf,
    m_sur = validation_data$surgery_mitosis,
    predicted_median = validation_data$predicted_median,
    location = validation_data$location,
    tumor_size = validation_data$size_mm
  )
  
  # order by surface
  plot_data_surface <- plot_data[order(plot_data$surface), ]
  plot_data_surface$x_position <- 1:nrow(plot_data_surface)
  
  # location colors
  location_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
  
  # plot
  jpeg("validation_output1/figures/figure2a_validation.jpg", 
       units = "in", width = 16, height = 10, res = 400)
  
  par(mfrow = c(2, 1), mar = c(10, 8, 6, 4), xaxs = "r", yaxs = "r")
  
  ### PANEL A: Signed difference ###
  plot(plot_data_surface$x_position, plot_data_surface$difference, 
       type = "h",
       lwd = 2.5,
       col = ifelse(plot_data_surface$difference > 0, 
                    alpha("#2ECC71", 0.6), 
                    alpha("#E74C3C", 0.6)),
       xlab = "", 
       ylab = "Signed Difference (Predicted - Observed)",
       main = "A) Prediction Error by Biopsy Surface Area",
       bty = "l",
       ylim = range(plot_data_surface$difference) * 1.15,
       xaxt = "n",
       cex.lab = 1.4,
       cex.main = 1.6,
       cex.axis = 1.2)
  
  abline(h = 0, lty = 1, col = "black", lwd = 2)
  
  # surface as line
  par(new = TRUE)
  plot(plot_data_surface$x_position, plot_data_surface$surface,
       type = "l",
       lwd = 3,
       col = alpha("navyblue", 0.7),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "",
       ylim = c(0, max(plot_data_surface$surface) * 1.1),
       bty = "n")
  
  axis(4, col = "navyblue", col.axis = "navyblue", cex.axis = 1.2)
  mtext("Biopsy Surface (HPF)", side = 4, line = 3, cex = 1.2, col = "navyblue")
  
  axis(1, at = plot_data_surface$x_position, 
       labels = plot_data_surface$patient_id,
       las = 2, cex.axis = 0.6)
  mtext("Patient ID", side = 1, line = 8, cex = 1.2)
  
  legend("topleft",
         legend = c("Overestimation", "Underestimation", "Biopsy surface"),
         col = c("#2ECC71", "#E74C3C", "navyblue"),
         lwd = c(3, 3, 3),
         lty = c(1, 1, 1),
         cex = 1.2, bty = "n")
  
  ### PANEL B: Absolute difference ###
  plot(plot_data_surface$x_position, plot_data_surface$abs_difference, 
       type = "h",
       lwd = 2.5,
       col = alpha("steelblue", 0.6),
       xlab = "", 
       ylab = "Absolute Difference |Predicted - Observed|",
       main = "B) Magnitude of Prediction Error",
       bty = "l",
       ylim = c(0, max(plot_data_surface$abs_difference) * 1.15),
       xaxt = "n",
       cex.lab = 1.4,
       cex.main = 1.6,
       cex.axis = 1.2)
  
  # Mean absolute error line
  abline(h = mean(plot_data_surface$abs_difference), 
         col = "red", lwd = 2.5, lty = 2)
  
  # Add surface as line
  par(new = TRUE)
  plot(plot_data_surface$x_position, plot_data_surface$surface,
       type = "l",
       lwd = 3,
       col = alpha("navyblue", 0.7),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "",
       ylim = c(0, max(plot_data_surface$surface) * 1.1),
       bty = "n")
  
  axis(4, col = "navyblue", col.axis = "navyblue", cex.axis = 1.2)
  mtext("Biopsy Surface (HPF)", side = 4, line = 3, cex = 1.2, col = "navyblue")
  
  axis(1, at = plot_data_surface$x_position, 
       labels = plot_data_surface$patient_id,
       las = 2, cex.axis = 0.6)
  mtext("Patient ID", side = 1, line = 8, cex = 1.2)
  
  legend("topleft",
         legend = c(sprintf("MAE = %.2f mitoses", 
                            mean(plot_data_surface$abs_difference)),
                    "Biopsy surface"),
         col = c("red", "navyblue"),
         lwd = c(2.5, 3),
         lty = c(2, 1),
         cex = 1.2, bty = "n")
  
  dev.off()
  
  cat("✓ Figure 2A saved: validation_output1/figures/figure2a_validation.jpg\n")
}

################################################################################
# FIGURE 2B: Prediction Error by Tumor Size
################################################################################

create_figure2b_validation <- function(validation_data, predicted_samples) {
  
  cat("\nCreating Figure 2B: Prediction Error by Tumor Size...\n")
  
  # calculate differences
  difference <- validation_data$predicted_median - validation_data$surgery_mitosis
  abs_difference <- abs(difference)
  
  # create data frame
  plot_data <- data.frame(
    patient_id = 1:nrow(validation_data),
    difference = difference,
    abs_difference = abs_difference,
    tumor_size = validation_data$size_mm,
    location = validation_data$location
  )
  
  # order by tumor size
  plot_data_ordered <- plot_data %>% arrange(tumor_size)
  plot_data_ordered$x_position <- 1:nrow(plot_data_ordered)
  
  # location labels and colors
  location_names <- c("Esophagus", "Stomach", "Duodenum", 
                      "Small intestine", "Colon-rectum")
  location_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
  
  # plot
  jpeg("validation_output1/figures/figure2b_validation.jpg", 
       units = "in", width = 18, height = 10, res = 400)
  
  par(mar = c(8, 8, 4, 4))
  
  plot(plot_data_ordered$x_position, plot_data_ordered$difference, 
       type = "n",
       xlab = "", 
       ylab = "Difference (Predicted - Observed Mitotic Count)",
       main = "Prediction Error by Tumor Size - External Validation",
       bty = "l",
       ylim = range(plot_data_ordered$difference) * 1.15,
       xaxt = "n",
       cex.lab = 1.5,
       cex.main = 1.7,
       cex.axis = 1.3)
  
  # background bars by location
  for (i in 1:nrow(plot_data_ordered)) {
    rect(xleft = plot_data_ordered$x_position[i] - 0.5,
         xright = plot_data_ordered$x_position[i] + 0.5,
         ybottom = par("usr")[3],
         ytop = par("usr")[4],
         col = alpha(location_colors[plot_data_ordered$location[i]], 0.2),
         border = NA)
  }
  
  # horizontal grid
  abline(h = seq(floor(min(plot_data_ordered$difference)), 
                 ceiling(max(plot_data_ordered$difference)), by = 2),
         col = "gray85", lty = 1)
  
  # zero line
  abline(h = 0, lty = 1, col = "black", lwd = 2.5)
  
  # difference segments
  segments(x0 = plot_data_ordered$x_position, 
           y0 = 0, 
           y1 = plot_data_ordered$difference,
           col = ifelse(plot_data_ordered$difference > 0, 
                        "#2ECC71", 
                        "#E74C3C"),
           lwd = 2.5)
  
  # X-axis: Tumor sizes
  axis(1, at = plot_data_ordered$x_position, 
       labels = round(plot_data_ordered$tumor_size, 0),
       las = 2,
       cex.axis = 0.6,
       col.axis = "darkgreen",
       col.ticks = "darkgreen")
  mtext("Tumor Size (mm)", side = 1, line = 4, cex = 1.3, col = "darkgreen")
  
  # legends
  legend("topright", 
         title = expression(bold("Prediction Direction")),
         legend = c("Overestimation (λ > y)", 
                    "Underestimation (λ < y)", 
                    "Perfect prediction"),
         col = c("#2ECC71", "#E74C3C", "black"),
         lwd = 3, cex = 1.2, bty = "n")
  
  legend("topleft",
         title = expression(bold("Tumor Location")),
         legend = location_names[unique(plot_data_ordered$location)],
         fill = alpha(location_colors[unique(plot_data_ordered$location)], 0.2),
         border = "black",
         cex = 1.1,
         bty = "n")
  
  dev.off()
  
  cat("✓ Figure 2B saved: validation_output1/figures/figure2b_validation.jpg\n")
}

################################################################################
# FIGURE 3: Risk Classification Performance
################################################################################

create_figure3_validation <- function(validation_data) {
  
  cat("\nCreating Figure 3: Risk Classification...\n")
  
  # load risk function
  source("risk_functions.R")
  
  # calculate risk categories
  risk_actual <- mapply(risk_str,
                        size = validation_data$size_mm,
                        mic = validation_data$surgery_mitosis,
                        site = validation_data$location)
  
  risk_predicted <- mapply(risk_str,
                           size = validation_data$size_mm,
                           mic = round(validation_data$predicted_median),
                           site = validation_data$location)
  
  risk_biopsy <- mapply(risk_str,
                        size = validation_data$size_mm,
                        mic = validation_data$biopsy_mitosis,
                        site = validation_data$location)
  
  # risk levels
  risk_levels <- c("none", "very_low", "low", "moderate", "high")
  
  risk_actual <- factor(risk_actual, levels = risk_levels)
  risk_predicted <- factor(risk_predicted, levels = risk_levels)
  risk_biopsy <- factor(risk_biopsy, levels = risk_levels)
  
  # calculate accuracies
  acc_model <- sum(risk_predicted == risk_actual) / length(risk_actual)
  acc_biopsy <- sum(risk_biopsy == risk_actual) / length(risk_actual)
  
  # plot
  jpeg("validation_output1/figures/figure3_validation.jpg", 
       units = "in", width = 12, height = 6, res = 300)
  
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
  
  ### Panel A: Accuracy comparison ###
  barplot(c(acc_biopsy, acc_model) * 100,
          names.arg = c("Biopsy\nAlone", "PROMETheus\nPrediction"),
          ylim = c(0, 100),
          ylab = "Classification Accuracy (%)",
          main = "A) Risk Classification Performance",
          col = c("coral", "seagreen3"),
          border = NA,
          cex.names = 1.2,
          cex.axis = 1.2,
          cex.lab = 1.3,
          cex.main = 1.4)
  
  text(c(0.7, 1.9), 
       c(acc_biopsy, acc_model) * 100 + 3,
       labels = sprintf("%.1f%%", c(acc_biopsy, acc_model) * 100),
       font = 2, cex = 1.3)
  
  mtext(sprintf("N = %d patients", nrow(validation_data)), 
        side = 3, line = 0.5, cex = 0.9)
  
  ### Panel B: Stratified by actual risk ###
  acc_by_risk <- data.frame(
    risk = risk_levels,
    n = as.numeric(table(risk_actual)),
    biopsy = sapply(risk_levels, function(r) {
      idx <- risk_actual == r
      if(sum(idx) > 0) sum(risk_biopsy[idx] == risk_actual[idx]) / sum(idx) else NA
    }),
    model = sapply(risk_levels, function(r) {
      idx <- risk_actual == r
      if(sum(idx) > 0) sum(risk_predicted[idx] == risk_actual[idx]) / sum(idx) else NA
    })
  )
  
  # remove empty categories
  acc_by_risk <- acc_by_risk[acc_by_risk$n > 0, ]
  
  acc_matrix <- t(as.matrix(acc_by_risk[, c("biopsy", "model")]))
  colnames(acc_matrix) <- acc_by_risk$risk
  
  bp <- barplot(acc_matrix * 100,
                beside = TRUE,
                ylim = c(0, 110),
                ylab = "Accuracy (%)",
                main = "B) Performance by Actual Risk Group",
                col = c("coral", "seagreen3"),
                border = NA,
                las = 2,
                cex.axis = 1.2,
                cex.lab = 1.3,
                cex.main = 1.4,
                legend.text = c("Biopsy", "Model"),
                args.legend = list(x = "topright", bty = "n", cex = 1.1))
  
  # Add sample sizes
  for(j in 1:ncol(acc_matrix)) {
    text(mean(bp[, j]), -8,
         labels = sprintf("n=%d", acc_by_risk$n[j]),
         cex = 0.9, xpd = TRUE)
  }
  
  dev.off()
  
  cat("✓ Figure 3 saved: validation_output1/figures/figure3_validation.jpg\n")
}

################################################################################
# CREATE ALL VALIDATION FIGURES
################################################################################

create_all_validation_figures <- function() {
  
  cat("CREATING VALIDATION FIGURES FOR THESIS")
 
  # Load results if not in environment
  if(!exists("results")) {
    results <- readRDS("validation_output1/results/validation_results.rds")
  }
  
  validation_data <- results$data
  predicted_samples <- results$posterior_samples
  
  # Create all figures
  create_figure1_validation(validation_data, predicted_samples)
  create_figure2a_validation(validation_data, predicted_samples)
  create_figure2b_validation(validation_data, predicted_samples)
  create_figure3_validation(validation_data)
  
  
  cat("ALL VALIDATION FIGURES CREATED")
  cat("Figures saved to:")
  cat("validation_output1/figures/")
}


# after running validation:
create_all_validation_figures()
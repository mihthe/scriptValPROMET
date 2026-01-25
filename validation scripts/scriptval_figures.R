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
  
  N <- nrow(plot_data_surface)
  
  # plot
  jpeg("validation_output1/figures/figure2a_validation.jpg", 
       units = "in", width = 16, height = 12, res = 400)
  
  # separate plot and legend areas
  layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1), heights = c(4, 4, 1.2))
  
  ### PANEL A: Signed difference ###
  par(mar = c(5, 9, 4, 6))
  
  plot(plot_data_surface$x_position, plot_data_surface$difference, 
       type = "h",
       lwd = 2.5,
       col = ifelse(plot_data_surface$difference > 0, 
                    alpha("#2ECC71", 0.6), 
                    alpha("#E74C3C", 0.6)),
       xlab = "", 
       ylab = "",  # add this manually with mtext for better control
       main = "A) Prediction Error by Biopsy Surface Area",
       bty = "l",
       ylim = range(plot_data_surface$difference) * 1.15,
       xaxt = "n",
       cex.lab = 1.4,
       cex.main = 1.6,
       cex.axis = 1.2)
  
  # add y-axis label with mtext for precise position
  mtext("Signed Difference\n(Predicted - Observed)", 
        side = 2, line = 5, cex = 1.2)
  
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
  mtext("Biopsy Surface (HPF)", side = 4, line = 3.5, cex = 1.1, col = "navyblue")
  
  # x-axis shows position (1 to N), with patient ID as secondary info
  axis(1, at = plot_data_surface$x_position, 
       labels = plot_data_surface$x_position,  # show position 1, 2, 3... (ordered by surface)
       las = 1, cex.axis = 0.7)
  mtext("Patients (ordered by increasing biopsy surface)", side = 1, line = 3, cex = 1.1)
  
  
  ### PANEL B: Absolute difference ###
  par(mar = c(5, 9, 4, 6))
  
  plot(plot_data_surface$x_position, plot_data_surface$abs_difference, 
       type = "h",
       lwd = 2.5,
       col = alpha("steelblue", 0.6),
       xlab = "", 
       ylab = "",
       main = "B) Magnitude of Prediction Error",
       bty = "l",
       ylim = c(0, max(plot_data_surface$abs_difference) * 1.15),
       xaxt = "n",
       cex.lab = 1.4,
       cex.main = 1.6,
       cex.axis = 1.2)
  
  mtext("Absolute Difference\n|Predicted - Observed|", 
        side = 2, line = 5, cex = 1.2)
  
  # mean absolute error line
  mae_value <- mean(plot_data_surface$abs_difference)
  abline(h = mae_value, col = "red", lwd = 2.5, lty = 2)
  
  # add surface as line
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
  mtext("Biopsy Surface (HPF)", side = 4, line = 3.5, cex = 1.1, col = "navyblue")
  
  axis(1, at = plot_data_surface$x_position, 
       labels = plot_data_surface$x_position,
       las = 1, cex.axis = 0.7)
  mtext("Patients (ordered by increasing biopsy surface)", side = 1, line = 3, cex = 1.1)
  
  
  ### PANEL C: LEGEND AREA (separate from plots) ###
  par(mar = c(0, 2, 0, 2))
  
  plot.new()
  
  # legend for Panel A (sx)
  legend("left",
         title = expression(bold("Panel A - Direction")),
         legend = c("Overestimation (model > actual)", 
                    "Underestimation (model < actual)", 
                    "Biopsy surface"),
         col = c("#2ECC71", "#E74C3C", "navyblue"),
         lwd = c(3, 3, 3),
         lty = c(1, 1, 1),
         cex = 1.1, 
         bty = "n",
         horiz = FALSE)
  
  # legend for Panel B (dx)
  legend("right",
         title = expression(bold("Panel B - Magnitude")),
         legend = c(sprintf("MAE = %.2f mitoses", mae_value),
                    "Absolute error",
                    "Biopsy surface"),
         col = c("red", "steelblue", "navyblue"),
         lwd = c(2.5, 3, 3),
         lty = c(2, 1, 1),
         cex = 1.1, 
         bty = "n",
         horiz = FALSE)
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
  
  # separate regions for plot and legend
  # matrix: row 1 = plot (height 4), row 2 = legend (height 1)
  layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(4, 1))
  
  # PANEL 1: Main plot
  par(mar = c(6, 6, 4, 2))
  
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
  
  # PANEL 2: Legend area
  par(mar = c(1, 6, 1, 2))
  
  # empty plot for legend
  plot.new()
  
  # legend 1: Prediction Direction (left side)
  legend("left",
         title = expression(bold("Prediction Direction")),
         legend = c("Overestimation (λ > y)", 
                    "Underestimation (λ < y)", 
                    "Perfect prediction"),
         col = c("#2ECC71", "#E74C3C", "black"),
         lwd = 3, 
         cex = 1.2, 
         bty = "n",
         horiz = FALSE)
  
  # legend 2: Tumor Location (right side)
  legend("right",
         title = expression(bold("Tumor Location")),
         legend = location_names[unique(plot_data_ordered$location)],
         fill = alpha(location_colors[unique(plot_data_ordered$location)], 0.2),
         border = "black",
         cex = 1.1,
         bty = "n",
         horiz = FALSE)
  
  dev.off()
  
  cat("✓ Figure 2B (v2) saved: validation_output1/figures/figure2b_validation_v2.jpg\n")
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
  
  # keep as character vectors for comparison
  risk_actual_chr <- as.character(risk_actual)
  risk_predicted_chr <- as.character(risk_predicted)
  risk_biopsy_chr <- as.character(risk_biopsy)
  
  # remove if NA cases
  valid_idx <- !is.na(risk_actual_chr) & !is.na(risk_predicted_chr) & !is.na(risk_biopsy_chr)
  
  if(sum(!valid_idx) > 0) {
    cat("Warning:", sum(!valid_idx), "cases excluded due to NA risk values\n")
  }
  
  risk_actual_valid <- risk_actual_chr[valid_idx]
  risk_predicted_valid <- risk_predicted_chr[valid_idx]
  risk_biopsy_valid <- risk_biopsy_chr[valid_idx]
  
  n_valid <- length(risk_actual_valid)
  
  # calculate accuracies
  matches_model <- risk_predicted_valid == risk_actual_valid
  matches_biopsy <- risk_biopsy_valid == risk_actual_valid
  
  acc_model <- sum(matches_model) / n_valid
  acc_biopsy <- sum(matches_biopsy) / n_valid
   
  # plot
  jpeg("validation_output1/figures/figure3_validation.jpg", 
       units = "in", width = 14, height = 7, res = 300)
  
  par(mfrow = c(1, 2), mar = c(6, 6, 4, 2))
  
  ### Panel A: Accuracy comparison ###
  
  # set bar heights and check they're not NA
  bar_heights_A <- c(acc_biopsy * 100, acc_model * 100)
  
  # replace any NA with 0 for plotting
  bar_heights_A[is.na(bar_heights_A)] <- 0
  
  par(mgp = c(3, 1.8, 0))
  
  bp_A <- barplot(bar_heights_A,
                  names.arg = c("Biopsy\nAlone", "PROMETheus\nPrediction"),
                  ylim = c(0, 110),
                  ylab = "Classification Accuracy (%)",
                  main = "A) Risk Classification Performance",
                  col = c("coral", "seagreen3"),
                  border = NA,
                  cex.names = 1.2,
                  cex.axis = 1.2,
                  cex.lab = 1.3,
                  cex.main = 1.4)
  
  # percentage labels on top of bars
  text(bp_A, 
       bar_heights_A + 5,
       labels = sprintf("%.1f%%", bar_heights_A),
       font = 2, cex = 1.3)
  
  # sample size
  mtext(sprintf("N = %d patients", n_valid), 
        side = 3, line = 0.5, cex = 0.9)
  
  ### Panel B: Stratified by actual risk ###
  
  # risk levels in order
  risk_levels <- c("none", "very_low", "low", "moderate", "high")
  
  # accuracy for each risk group
  acc_by_risk <- data.frame(
    risk = risk_levels,
    n = sapply(risk_levels, function(r) sum(risk_actual_valid == r)),
    biopsy = sapply(risk_levels, function(r) {
      idx <- which(risk_actual_valid == r)
      if(length(idx) > 0) {
        sum(risk_biopsy_valid[idx] == risk_actual_valid[idx]) / length(idx) * 100
      } else {
        NA
      }
    }),
    model = sapply(risk_levels, function(r) {
      idx <- which(risk_actual_valid == r)
      if(length(idx) > 0) {
        sum(risk_predicted_valid[idx] == risk_actual_valid[idx]) / length(idx) * 100
      } else {
        NA
      }
    })
  )
  
  # remove empty categories
  acc_by_risk <- acc_by_risk[acc_by_risk$n > 0, ]
  
  if(nrow(acc_by_risk) == 0) {
    plot.new()
    text(0.5, 0.5, "No valid risk categories to display", cex = 1.5)
  } else {
  
    # matrix for grouped barplot
    acc_matrix <- rbind(
      Biopsy = acc_by_risk$biopsy,
      Model = acc_by_risk$model
    )
    colnames(acc_matrix) <- acc_by_risk$risk
    
    # NA to 0 for plotting
    acc_matrix[is.na(acc_matrix)] <- 0
    
    # combined labels with line break to avoid overlap
    combined_labels <- paste0(acc_by_risk$risk, "\n(n=", acc_by_risk$n, ")")
    
    #par(mar = c(8, 5, 4, 2))
    par(mar = c(7, 5, 4, 2))
    
    bp_B <- barplot(acc_matrix,
                    beside = TRUE,
                    ylim = c(0, 120),
                    ylab = "Accuracy (%)",
                    main = "B) Performance by Actual Risk Group",
                    col = c("coral", "seagreen3"),
                    border = NA,
                    names.arg = combined_labels,  # use combined labels directly
                    cex.names = 1.0,              # adjust size if needed
                    cex.axis = 1.2,
                    cex.lab = 1.3,
                    cex.main = 1.4,
                    legend.text = c("Biopsy", "Model"),
                    args.legend = list(x = "topright", bty = "n", cex = 1.1))
  }
    
  #   # bp_B returns a matrix: rows = groups (Biopsy, Model), cols = risk categories
  #   # center position for each risk category
  #   category_centers <- colMeans(bp_B)
  #   
  #   # risk category names (line 1 below axis)
  #   axis(1, at = category_centers, 
  #        labels = acc_by_risk$risk,
  #        tick = FALSE,
  #        line = 0.5,
  #        cex.axis = 1.1)
  #   
  #   # sample sizes below the category names (line 2)
  #   mtext(sprintf("(n=%d)", acc_by_risk$n),
  #         side = 1,
  #         at = category_centers,
  #         line = 2.5,
  #         cex = 0.9)
  #   
  #   # add x-axis title
  #   mtext("Actual Risk Category", side = 1, line = 5, cex = 1.2)
  # }
  
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

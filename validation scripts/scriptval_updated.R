## updated validation script

# Libraries
library(rethinking)
library(ggplot2)
library(dplyr)
library(caret)

################################################################################
# STEP 1: LOAD TRAINED MODEL AND STANDARDIZATION PARAMETERS
################################################################################

# run data simulation file to
# get the trained model (m) 
# and posterior samples (post)
# and save 
  # posterior samples
  saveRDS(post, "validation_output1/model_posterior.rds")
  
  # standardization parameters from training data
  std_params <- list(
    Si_mean = mean(sim$Si),
    Si_sd = sd(sim$Si),
    m_bio_mean = attr(dat$m_bio, 'scaled:center'),
    m_bio_sd = attr(dat$m_bio, 'scaled:scale'),
    Su_mean = attr(dat$Su, 'scaled:center'),
    Su_sd = attr(dat$Su, 'scaled:scale')
  )

  saveRDS(std_params, "validation_output1/standardization_params.rds")
  
 

################################################################################
# STEP 2: LOAD AND PREPARE REAL VALIDATION DATA
################################################################################

load_validation_data <- function(filepath = "validationset.csv") {
  
  # read validation data
  validation_raw <- read.csv(filepath, stringsAsFactors = FALSE)
  
  # check required columns
  required_cols <- c("patient_id", "size_mm", "location", 
                     "biopsy_mitosis", "biopsy_surface_hpf", "surgery_mitosis")
  
  missing_cols <- setdiff(required_cols, names(validation_raw))
  if(length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", 
                 paste(missing_cols, collapse = ", ")))
  }
  
  # load standardization parameters from TRAINING data
  std_params <- readRDS("validation_output1/standardization_params.rds")
  
  # standardize using TRAINING data statistics 
  validation_data <- validation_raw %>%
    mutate(
      # standardize size (convert mm to training scale)
      size_raw = (size_mm - std_params$Si_mean) / std_params$Si_sd,
      
      # keep raw values for interpretation
      biopsy_mitosis_raw = biopsy_mitosis,
      biopsy_surface_raw = biopsy_surface_hpf,
      
      # add data source identifier
      data_source = "real_validation"
    )
  
  return(list(
    data = validation_data,
    std_params = std_params
  ))
}

################################################################################
# STEP 3: PREDICTION FUNCTION (USING TRAINED MODEL)
################################################################################

predict_mitosis <- function(patient_data, posterior, std_params) { 
  
  # standardize inputs using TRAINING parameters
  m_bio_std <- (patient_data$biopsy_mitosis - std_params$m_bio_mean) / 
    std_params$m_bio_sd
  Su_std <- (patient_data$biopsy_surface_hpf - std_params$Su_mean) / 
    std_params$Su_sd
  Si_std <- patient_data$size_raw  
  
  n_samples <- length(posterior$a)
  
  # compute log(lambda) for each posterior sample
  log_lambda <- numeric(n_samples)
  
  for(i in 1:n_samples) {
    log_lambda[i] <- posterior$a[i] + 
      posterior$b[i, patient_data$location] * Si_std +
      posterior$g[i] * m_bio_std +
      posterior$e[i, patient_data$location] * Su_std
  }
  
  # convert to expected counts
  lambda <- exp(log_lambda)
  
  # generate predicted counts
  predicted_counts <- rpois(n_samples, lambda = lambda)
  
  return(predicted_counts)
}

################################################################################
# STEP 4: PREDICT FOR ALL VALIDATION PATIENTS
################################################################################

predict_all_validation <- function(validation_data, posterior, std_params) { 
  
  n_patients <- nrow(validation_data)
  n_samples <- length(posterior$a)
  
  cat("MAKING PREDICTIONS\n")
  
  # storage for predictions
  predicted_samples <- matrix(NA, nrow = n_samples, ncol = n_patients)
  
  # predict for each patient
  for(i in 1:n_patients) {
    
    predicted_samples[, i] <- predict_mitosis(
      patient_data = validation_data[i, ],
      posterior = posterior,
      std_params = std_params
    )
    
    if(i %% 10 == 0 | i == n_patients) {
      cat(sprintf("  Progress: %d/%d\r", i, n_patients))
    }
  }
  
  cat("\n✓ Predictions complete\n\n")
  
  # calculate summary statistics
  validation_data <- validation_data %>%
    mutate(
      predicted_mean = colMeans(predicted_samples),
      predicted_median = apply(predicted_samples, 2, median),
      predicted_ci_lower = apply(predicted_samples, 2, function(x) HPDI(x, 0.95)[1]),
      predicted_ci_upper = apply(predicted_samples, 2, function(x) HPDI(x, 0.95)[2]),
      prediction_error = predicted_mean - surgery_mitosis,
      abs_error = abs(prediction_error),
      
      # calculate coverage
      in_ci = (surgery_mitosis >= predicted_ci_lower) & 
        (surgery_mitosis <= predicted_ci_upper)
    )
  
  return(list(
    data = validation_data,
    posterior_samples = predicted_samples
  ))
}

################################################################################
# STEP 5: CALCULATE VALIDATION METRICS
################################################################################

# Miettinen & Lasota Risk Classification
classify_risk_ML <- function(size_mm, mitotic_count, location) {
  size_cm <- size_mm / 10
  is_gastric <- (location == 2)
  
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

calculate_validation_metrics <- function(validation_results) {
  
  data <- validation_results$data
  
  cat("CALCULATING METRICS\n")

  # 1. continuous prediction metrics
  mae <- mean(abs(data$prediction_error))
  rmse <- sqrt(mean(data$prediction_error^2))
  bias <- mean(data$prediction_error)
  r_squared <- cor(data$predicted_mean, data$surgery_mitosis)^2
  
  # 2. coverage
  coverage_95 <- mean(data$in_ci)
  
  # 3. risk classification
  risk_levels <- c("none", "very_low", "low", "moderate", "high")
  
  # actual risk from surgery
  actual_risk <- mapply(classify_risk_ML,
                        size_mm = data$size_mm,
                        mitotic_count = data$surgery_mitosis,
                        location = data$location)
  
  # predicted risk from model
  predicted_risk <- mapply(classify_risk_ML,
                           size_mm = data$size_mm,
                           mitotic_count = round(data$predicted_mean),
                           location = data$location)
  
  # biopsy risk
  biopsy_risk <- mapply(classify_risk_ML,
                        size_mm = data$size_mm,
                        mitotic_count = data$biopsy_mitosis,
                        location = data$location)
  
  # convert to factors 
  actual_risk <- factor(actual_risk, levels = risk_levels)
  predicted_risk <- factor(predicted_risk, levels = risk_levels)
  biopsy_risk <- factor(biopsy_risk, levels = risk_levels)
  
  # confusion matrices
  cm_model <- confusionMatrix(predicted_risk, actual_risk)
  cm_biopsy <- confusionMatrix(biopsy_risk, actual_risk)
  
  # compile results
  metrics <- list(
    continuous = list(
      mae = mae,
      rmse = rmse,
      bias = bias,
      r_squared = r_squared,
      coverage_95 = coverage_95,
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
  
  # print summary
  cat(sprintf("MAE: %.2f mitoses\n", mae))
  cat(sprintf("RMSE: %.2f mitoses\n", rmse))
  cat(sprintf("Bias: %.2f mitoses\n", bias))
  cat(sprintf("R²: %.3f\n", r_squared))
  cat(sprintf("95%% CI Coverage: %.1f%%\n\n", coverage_95 * 100))
  
  cat(sprintf("Model Classification Accuracy: %.1f%%\n", 
              cm_model$overall["Accuracy"] * 100))
  cat(sprintf("Biopsy Classification Accuracy: %.1f%%\n", 
              cm_biopsy$overall["Accuracy"] * 100))
  
  return(metrics)
}

################################################################################
# STEP 6: CREATE VALIDATION PLOTS
################################################################################

create_validation_plots <- function(validation_results, metrics) {
  
  data <- validation_results$data
  
  cat("CREATING PLOTS\n")

  # Plot 1: Calibration plot
  pdf("validation_output1/figures/calibration_plot.pdf", width = 10, height = 8)
  
  par(mar = c(5, 5, 4, 2))
  
  max_val <- max(c(data$predicted_mean, data$surgery_mitosis))
  
  plot(data$predicted_mean, data$surgery_mitosis,
       xlab = "Predicted Mitotic Count (PROMETheus)",
       ylab = "Observed Mitotic Count (Surgery)",
       main = "External Validation: Calibration Plot",
       pch = 16, col = alpha(data$location, 0.6), cex = 1.5,
       xlim = c(0, max_val * 1.1), ylim = c(0, max_val * 1.1))
  
  # prediction line
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  
  # 95% CI bars
  segments(x0 = data$predicted_ci_lower,
           y0 = data$surgery_mitosis,
           x1 = data$predicted_ci_upper,
           y1 = data$surgery_mitosis,
           col = rgb(0, 0, 0, 0.2))
  
  legend("topleft",
         legend = c(
           sprintf("MAE = %.2f", metrics$continuous$mae),
           sprintf("RMSE = %.2f", metrics$continuous$rmse),
           sprintf("R² = %.3f", metrics$continuous$r_squared),
           sprintf("Coverage = %.1f%%", metrics$continuous$coverage_95 * 100),
           "Perfect prediction"
         ),
         lty = c(0, 0, 0, 0, 2),
         col = c("black", "black", "black", "black", "red"),
         lwd = c(0, 0, 0, 0, 2),
         bty = "n", cex = 1.1)
  
  dev.off()
  
  cat("✓ Calibration plot created\n")
  
  # Plot 2: Residuals
  pdf("validation_output1/figures/residuals_plot.pdf", width = 10, height = 8)
  
  par(mar = c(5, 5, 4, 2))
  
  plot(data$predicted_mean, data$prediction_error,
       xlab = "Predicted Mitotic Count",
       ylab = "Residual (Predicted - Observed)",
       main = "Residual Plot",
       pch = 16, col = alpha(data$location, 0.6), cex = 1.5)
  
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  abline(h = c(-2, 2), col = "orange", lwd = 1, lty = 3)
  
  dev.off()
  
  cat("✓ Residual plot created\n")
  
  # Plot 3: Risk classification comparison
  pdf("validation_output1/figures/risk_classification.pdf", width = 12, height = 6)
  
  par(mfrow = c(1, 2), mar = c(8, 5, 4, 2))
  
  # model confusion matrix
  cm_model_prop <- prop.table(metrics$classification$cm_model, margin = 1) * 100
  
  image(1:ncol(cm_model_prop), 1:nrow(cm_model_prop), 
        t(cm_model_prop),
        col = colorRampPalette(c("white", "steelblue"))(20),
        xlab = "Observed Risk", ylab = "Predicted Risk (Model)",
        main = "PROMETheus Classification",
        axes = FALSE)
  
  axis(1, at = 1:ncol(cm_model_prop), labels = colnames(cm_model_prop), las = 2)
  axis(2, at = 1:nrow(cm_model_prop), labels = rownames(cm_model_prop), las = 1)
  
  # add text
  for(i in 1:nrow(metrics$classification$cm_model)) {
    for(j in 1:ncol(metrics$classification$cm_model)) {
      text(j, i, 
           sprintf("%d\n(%.0f%%)", 
                   metrics$classification$cm_model[i,j], 
                   cm_model_prop[i,j]),
           col = ifelse(cm_model_prop[i,j] > 50, "white", "black"),
           cex = 1.2, font = 2)
    }
  }
  
  # biopsy confusion matrix
  cm_biopsy_prop <- prop.table(metrics$classification$cm_biopsy, margin = 1) * 100
  
  image(1:ncol(cm_biopsy_prop), 1:nrow(cm_biopsy_prop), 
        t(cm_biopsy_prop),
        col = colorRampPalette(c("white", "coral"))(20),
        xlab = "Observed Risk", ylab = "Predicted Risk (Biopsy)",
        main = "Biopsy Alone Classification",
        axes = FALSE)
  
  axis(1, at = 1:ncol(cm_biopsy_prop), labels = colnames(cm_biopsy_prop), las = 2)
  axis(2, at = 1:nrow(cm_biopsy_prop), labels = rownames(cm_biopsy_prop), las = 1)
  
  for(i in 1:nrow(metrics$classification$cm_biopsy)) {
    for(j in 1:ncol(metrics$classification$cm_biopsy)) {
      text(j, i, 
           sprintf("%d\n(%.0f%%)", 
                   metrics$classification$cm_biopsy[i,j], 
                   cm_biopsy_prop[i,j]),
           col = ifelse(cm_biopsy_prop[i,j] > 50, "white", "black"),
           cex = 1.2, font = 2)
    }
  }
  
  dev.off()
  
  cat("✓ Risk classification plots created\n\n")
}

################################################################################
# VALIDATION WORKFLOW
################################################################################

run_validation <- function() {
  
  # load trained model
  if(!file.exists("validation_output1/model_posterior.rds")) {
    stop("Model not found. Run data analysis file first and save the model.")
  }
  
  posterior <- readRDS("validation_output1/model_posterior.rds")
  cat("\n✓ Loaded trained model\n")
  
  # load validation data
  loaded <- load_validation_data("validationset.csv")
  validation_data <- loaded$data
  std_params <- loaded$std_params
  
  # make predictions
  validation_results <- predict_all_validation(validation_data, posterior, std_params)
  
  # calculate metrics
  metrics <- calculate_validation_metrics(validation_results)
  
  # create plots
  create_validation_plots(validation_results, metrics)
  
  # save results
  results <- list(
    data = validation_results$data,
    posterior_samples = validation_results$posterior_samples,
    metrics = metrics,
    info = list(
      n_patients = nrow(validation_data),
      date = Sys.Date()
    )
  )
  
  saveRDS(results, "validation_output1/results/validation_results.rds")
  write.csv(validation_results$data, 
            "validation_output1/results/predictions.csv", 
            row.names = FALSE)
  
  
  cat("VALIDATION COMPLETE")
  
  
  return(results)
}


# run validation
results <- run_validation()
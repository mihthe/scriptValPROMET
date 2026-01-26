### risk group differences - validation data

# run simulation file
# save posterior and standardisation parametres
# run validation dataset file
# run validation file

# libraries
library(rethinking)
library(readxl)
library(dplyr)

# output directory if it doesn't exist
#dir.create("validation_output1/figures", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PART 1: load data from validation file 
# =============================================================================

# load all required data
load_validation_data <- function() {
  
  # load posterior
  post <- readRDS("validation_output1/model_posterior.rds")
  cat("  Posterior samples loaded:", length(post$a), "samples\n")
  
  # load standardization parameters
  std_params <- readRDS("validation_output1/standardization_params.rds")
  cat("  Standardization parameters loaded\n")
  
  # load validation dataset
  validation_raw <- read_excel("validation_dataset_ready.xlsx")
  validation_raw <- as.data.frame(validation_raw)
  cat("  Validation data loaded:", nrow(validation_raw), "patients\n")
  
  # rename columns 
  # Si (size), m_bio, m_sur (surgery mitosis), Su (surface), L (location)
  validation_data <- validation_raw %>%
    rename(
      Si = size_mm,
      m_sur = m_surg,
      Su = biopsy_surface_hpf
    )
  
  return(list(
    data = validation_data,
    post = post,
    std_params = std_params
  ))
}

# =============================================================================
# PART 2: risk classification functions
# =============================================================================

# Miettinen & Lasota Risk Classification (single case)
risk_str <- function(size, mic, site) {
  # size: tumor size in mm
  # mic: mitotic count per 5mm²
  # site: location code (4 = stomach/gastric)
  
  size <- size / 10  # convert mm to cm
  
  if (site == 4L) {  # gastric
    if (mic <= 5) {
      if (size <= 2) return("none")
      else if (size <= 5) return("very low")
      else if (size <= 10) return("low")
      else return("moderate")
    } else {  # mic > 5
      if (size <= 2) return("none")
      else if (size <= 5) return("moderate")
      else return("high")
    }
  } else {  # Non-gastric
    if (mic <= 5) {
      if (size <= 2) return("none")
      else if (size <= 5) return("low")
      else if (size <= 10) return("moderate")
      else return("high")
    } else {  # mic > 5
      return("high")
    }
  }
}

# vectorized version
risk_str_vectorized <- function(size, mic, site) {
  n <- max(length(size), length(mic), length(site))
  risk <- character(n)
  
  for (i in 1:n) {
    risk[i] <- risk_str(
      size = size[i],
      mic = mic[i],
      site = site[i]
    )
  }
  return(risk)
}

# risk classifications using VALIDATION data
compute_risk_classifications <- function(validation_data, post, std_params) {
  # validation_data: data frame with Si, m_bio, m_sur, Su, L columns
  # post: posterior samples from fitted model
  # std_params: standardization parameters from training data
  
  n_cases <- nrow(validation_data)
  
  # 1. Risk based on BIOPSY mitotic count
  risk_biopsy <- risk_str_vectorized(
    size = validation_data$Si,
    mic = validation_data$m_bio,
    site = validation_data$L
  )
  
  # 2. Risk based on ACTUAL SURGERY mitotic count
  risk_surgery_actual <- risk_str_vectorized(
    size = validation_data$Si,
    mic = validation_data$m_sur,
    site = validation_data$L
  )
  
  # 3. Risk based on PREDICTED surgery mitotic count
  # Link function adapted for validation data
  m_link <- function(Si, m_bio, Su, L) {
    # Standardize using TRAINING parameters
    Si_std <- (Si - std_params$Si_mean) / std_params$Si_sd
    m_bio_std <- (m_bio - std_params$m_bio_mean) / std_params$m_bio_sd
    Su_std <- (Su - std_params$Su_mean) / std_params$Su_sd
    
    # log(lambda) using posterior samples
    mu <- with(post, {
      a + b[, L] * Si_std + g * m_bio_std + e[, L] * Su_std
    })
    lambda <- exp(mu)
    return(lambda)
  }
  
  # predictions for each case
  risk_surgery_predicted <- character(n_cases)
  predicted_mitosis_median <- numeric(n_cases)
  
  for (i in 1:n_cases) {
    # Lambda distribution for this case
    lambda <- m_link(
      Si = validation_data$Si[i],
      m_bio = validation_data$m_bio[i],
      Su = validation_data$Su[i],
      L = validation_data$L[i]
    )
    
    # median of lambda as predicted mitotic count
    pred_mitosis <- median(lambda)
    predicted_mitosis_median[i] <- pred_mitosis
    
    # risk class from predicted mitosis
    risk_surgery_predicted[i] <- risk_str(
      size = validation_data$Si[i],
      mic = round(pred_mitosis),
      site = validation_data$L[i]
    )
    
    if (i %% 10 == 0 || i == n_cases) {
      cat("\r  Progress:", i, "/", n_cases, "cases")
    }
  }
  cat("\n")
  
  # results
  result <- data.frame(
    case_id = 1:n_cases,
    size = validation_data$Si,
    location = validation_data$L,
    m_bio = validation_data$m_bio,
    m_sur = validation_data$m_sur,
    surface_bio = validation_data$Su,
    risk_biopsy = risk_biopsy,
    risk_surgery_actual = risk_surgery_actual,
    risk_surgery_predicted = risk_surgery_predicted,
    predicted_mitosis = predicted_mitosis_median,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# =============================================================================
# PART 3: stratified accuracy functions
# =============================================================================

# stratified accuracy by risk group
calculate_stratified_accuracy <- function(risk_data) {
  
  risk_levels <- c("none", "very low", "low", "moderate", "high")
  results_list <- list()
  
  for (risk_level in risk_levels) {
    subset_data <- risk_data[risk_data$risk_surgery_actual == risk_level, ]
    
    if (nrow(subset_data) > 0) {
      # model accuracy
      correct_model <- sum(subset_data$risk_surgery_predicted == subset_data$risk_surgery_actual)
      accuracy_model <- correct_model / nrow(subset_data)
      
      # biopsy accuracy
      correct_biopsy <- sum(subset_data$risk_biopsy == subset_data$risk_surgery_actual)
      accuracy_biopsy <- correct_biopsy / nrow(subset_data)
      
      results_list[[risk_level]] <- data.frame(
        Risk_Group = risk_level,
        N_Cases = nrow(subset_data),
        Biopsy_Accuracy = accuracy_biopsy,
        Model_Accuracy = accuracy_model,
        Improvement = accuracy_model - accuracy_biopsy,
        stringsAsFactors = FALSE
      )
    }
  }
  
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  return(results)
}

# stratified accuracy by biopsy surface
calculate_surface_stratified_accuracy <- function(risk_data) {
  
  # surface categories
  risk_data$surface_category <- cut(
    risk_data$surface_bio,
    breaks = c(0, 5, 10, 15, 20, 23.5),
    labels = c("0-5", "5-10", "10-15", "15-20", "20-23.5"),
    include.lowest = TRUE
  )
  
  results_list <- list()
  surface_cats <- levels(risk_data$surface_category)
  
  for (cat in surface_cats) {
    subset_data <- risk_data[risk_data$surface_category == cat, ]
    
    if (nrow(subset_data) > 0) {
      # model accuracy
      correct_model <- sum(subset_data$risk_surgery_predicted == subset_data$risk_surgery_actual,
                           na.rm = TRUE)
      accuracy_model <- correct_model / nrow(subset_data)
      
      # biopsy accuracy
      correct_biopsy <- sum(subset_data$risk_biopsy == subset_data$risk_surgery_actual,
                            na.rm = TRUE)
      accuracy_biopsy <- correct_biopsy / nrow(subset_data)
      
      # mean values
      mean_predicted <- mean(subset_data$predicted_mitosis, na.rm = TRUE)
      mean_actual <- mean(subset_data$m_sur, na.rm = TRUE)
      
      results_list[[cat]] <- data.frame(
        Surface_Category = cat,
        N_Cases = nrow(subset_data),
        Mean_Surface = mean(subset_data$surface_bio, na.rm = TRUE),
        Biopsy_Accuracy = accuracy_biopsy,
        Model_Accuracy = accuracy_model,
        Improvement = accuracy_model - accuracy_biopsy,
        Mean_Predicted_Mitosis = mean_predicted,
        Mean_Actual_Mitosis = mean_actual,
        stringsAsFactors = FALSE
      )
    }
  }
  
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  return(results)
}

# =============================================================================
# PART 4: plots
# =============================================================================

# -----------------------------------------------------------------------------
# Plot 1: accuracy by risk group 
# -----------------------------------------------------------------------------
plot_accuracy_by_risk_group <- function(stratified_results, save_plot = TRUE,
                                        output_dir = "validation_output1/figures/") {
  
  if (save_plot) {
    jpeg(paste0(output_dir, "accuracy_by_risk_group.jpg"),
         units = "in", width = 12, height = 8, res = 300)
  }
  
  # margins 
  par(mfrow = c(1, 1), mar = c(10, 6, 4, 2))
  
  # data
  acc_matrix <- t(as.matrix(stratified_results[, c("Biopsy_Accuracy", "Model_Accuracy")]))
  colnames(acc_matrix) <- stratified_results$Risk_Group
  
  # barplot
  bp <- barplot(
    acc_matrix * 100,
    beside = TRUE,
    ylim = c(0, 120),  
    ylab = "",
    xlab = "",
    main = "Risk Classification Accuracy by Actual Risk Group\n(External Validation)",
    col = c("skyblue", "seagreen3"),
    border = NA,
    xaxt = "n",  
    cex.main = 1.3
  )
  
  # y-axis label 
  mtext("Accuracy (%)", side = 2, line = 4, cex = 1.2, font = 2)
  
  # x-axis labels at 45 degrees to avoid overlap
  x_positions <- colMeans(bp)
  text(x = x_positions, 
       y = par("usr")[3] - 5,  
       labels = stratified_results$Risk_Group,
       srt = 45,  
       adj = 1,   
       xpd = TRUE,
       cex = 1.1,
       font = 2)
  
  # x-axis title 
  mtext("Actual Risk Group (Surgery)", side = 1, line = 6.5, cex = 1.2, font = 2)
  
  # % labels on bars
  for (i in 1:nrow(acc_matrix)) {
    for (j in 1:ncol(acc_matrix)) {
      text(
        x = bp[i, j],
        y = acc_matrix[i, j] * 100 + 4,
        labels = paste0(round(acc_matrix[i, j] * 100, 1), "%"),
        font = 2,
        cex = 0.9
      )
    }
  }
  
  # sample sizes below the bar groups
  for (j in 1:ncol(acc_matrix)) {
    text(
      x = mean(bp[, j]),
      y = -12,
      labels = paste0("n=", stratified_results$N_Cases[j]),
      cex = 0.9,
      xpd = TRUE,
      font = 3
    )
  }
  
  # legend
  legend("topright",
         legend = c("Biopsy alone", "PROMETheus model"),
         fill = c("skyblue", "seagreen3"),
         border = NA,
         bty = "n",
         cex = 1.1)
  
  if (save_plot) dev.off()
  
}

# -----------------------------------------------------------------------------
# Plot 2: accuracy by biopsy surface 
# -----------------------------------------------------------------------------
plot_accuracy_by_surface <- function(surface_results, save_plot = TRUE,
                                     output_dir = "validation_output1/figures/") {
  
  if (save_plot) {
    jpeg(paste0(output_dir, "accuracy_by_surface.jpg"),
         units = "in", width = 12, height = 8, res = 300)
  }
  
  # margins
  par(mfrow = c(1, 1), mar = c(10, 6, 4, 2))
  
  # data
  acc_matrix <- t(as.matrix(surface_results[, c("Biopsy_Accuracy", "Model_Accuracy")]))
  colnames(acc_matrix) <- surface_results$Surface_Category
  
  # barplot
  bp <- barplot(
    acc_matrix * 100,
    beside = TRUE,
    ylim = c(0, 120),
    ylab = "",
    xlab = "",
    main = "Risk Classification Accuracy by Biopsy Surface Size\n(External Validation)",
    col = c("skyblue", "seagreen3"),
    border = NA,
    xaxt = "n",
    cex.main = 1.3
  )
  
  # y-axis label
  mtext("Accuracy (%)", side = 2, line = 4, cex = 1.2, font = 2)
  
  # x-axis labels 
  x_positions <- colMeans(bp)
  surface_labels <- paste0(surface_results$Surface_Category, " HPF")
  text(x = x_positions,
       y = par("usr")[3] - 5,
       labels = surface_labels,
       srt = 45,
       adj = 1,
       xpd = TRUE,
       cex = 1.0,
       font = 2)
  
  # x-axis title
  mtext("Biopsy Surface (High Power Fields)", side = 1, line = 6.5, cex = 1.2, font = 2)
  
  # % labels
  for (i in 1:nrow(acc_matrix)) {
    for (j in 1:ncol(acc_matrix)) {
      text(
        x = bp[i, j],
        y = acc_matrix[i, j] * 100 + 4,
        labels = paste0(round(acc_matrix[i, j] * 100, 1), "%"),
        font = 2,
        cex = 0.85
      )
    }
  }
  
  # sample sizes
  for (j in 1:ncol(acc_matrix)) {
    text(
      x = mean(bp[, j]),
      y = -12,
      labels = paste0("n=", surface_results$N_Cases[j]),
      cex = 0.9,
      xpd = TRUE,
      font = 3
    )
  }
  
  # legend
  legend("topright",
         legend = c("Biopsy alone", "PROMETheus model"),
         fill = c("skyblue", "seagreen3"),
         border = NA,
         bty = "n",
         cex = 1.1)
  
  if (save_plot) dev.off()
  
}

# -----------------------------------------------------------------------------
# Plot 3: improvement by surface 
# -----------------------------------------------------------------------------
plot_improvement_by_surface <- function(surface_results, save_plot = TRUE,
                                        output_dir = "validation_output1/figures/") {
  
  if (save_plot) {
    jpeg(paste0(output_dir, "improvement_by_surface.jpg"),
         units = "in", width = 10, height = 8, res = 300)
  }
  
  # margins
  par(mfrow = c(1, 1), mar = c(10, 6, 4, 2))
  
  # y-axis limits
  y_max <- max(surface_results$Improvement * 100, 10, na.rm = TRUE) + 15
  y_min <- min(surface_results$Improvement * 100, -10, na.rm = TRUE) - 10
  
  # barplot
  bp <- barplot(
    surface_results$Improvement * 100,
    ylim = c(y_min, y_max),
    ylab = "",
    xlab = "",
    main = "Model Improvement Over Biopsy Alone\nby Biopsy Surface Size (External Validation)",
    col = ifelse(surface_results$Improvement > 0, "seagreen3", "coral"),
    border = NA,
    xaxt = "n",
    cex.main = 1.3
  )
  
  # reference line at zero
  abline(h = 0, lty = 2, lwd = 2, col = "gray40")
  
  # y-axis label
  mtext("Accuracy Improvement (percentage points)", side = 2, line = 4, cex = 1.1, font = 2)
  
  # x-axis labels with rotation
  surface_labels <- paste0(surface_results$Surface_Category, " HPF")
  text(x = bp,
       y = par("usr")[3] - (y_max - y_min) * 0.05,
       labels = surface_labels,
       srt = 45,
       adj = 1,
       xpd = TRUE,
       cex = 1.0,
       font = 2)
  
  # x-axis title
  mtext("Biopsy Surface (High Power Fields)", side = 1, line = 6.5, cex = 1.2, font = 2)
  
  # value labels on bars
  label_offset <- ifelse(surface_results$Improvement > 0, 3, -3)
  text(
    x = bp,
    y = surface_results$Improvement * 100 + label_offset,
    labels = paste0(ifelse(surface_results$Improvement > 0, "+", ""),
                    round(surface_results$Improvement * 100, 1), "%"),
    font = 2,
    cex = 1.0
  )
  
  # sample sizes at bottom
  for (j in 1:length(bp)) {
    text(
      x = bp[j],
      y = y_min + (y_max - y_min) * 0.02,
      labels = paste0("n=", surface_results$N_Cases[j]),
      cex = 0.9,
      xpd = TRUE,
      font = 3
    )
  }
  
  # legend
  legend("topright",
         legend = c("Model better than biopsy", "Biopsy better than model"),
         fill = c("seagreen3", "coral"),
         border = NA,
         bty = "n",
         cex = 1.0)
  
  if (save_plot) dev.off()
  
  cat("  Plot saved: improvement_by_surface.jpg\n")
}

# =============================================================================
# PART 5: accuracy metrics overall
# =============================================================================

calculate_overall_accuracy <- function(risk_data) {
  
  n <- nrow(risk_data)
  
  # overall accuracy
  model_correct <- sum(risk_data$risk_surgery_predicted == risk_data$risk_surgery_actual)
  biopsy_correct <- sum(risk_data$risk_biopsy == risk_data$risk_surgery_actual)
  
  model_accuracy <- model_correct / n
  biopsy_accuracy <- biopsy_correct / n
  
  # high-risk specific metrics
  actual_high <- risk_data$risk_surgery_actual == "high"
  
  # model: sensitivity and specificity for high risk
  model_pred_high <- risk_data$risk_surgery_predicted == "high"
  model_sensitivity <- sum(actual_high & model_pred_high) / sum(actual_high)
  model_specificity <- sum(!actual_high & !model_pred_high) / sum(!actual_high)
  
  # biopsy: sensitivity and specificity for high risk
  biopsy_pred_high <- risk_data$risk_biopsy == "high"
  biopsy_sensitivity <- sum(actual_high & biopsy_pred_high) / sum(actual_high)
  biopsy_specificity <- sum(!actual_high & !biopsy_pred_high) / sum(!actual_high)
  
  results <- list(
    overall = data.frame(
      Method = c("Biopsy alone", "PROMETheus model"),
      Correct = c(biopsy_correct, model_correct),
      Total = c(n, n),
      Accuracy = c(biopsy_accuracy, model_accuracy),
      stringsAsFactors = FALSE
    ),
    high_risk = data.frame(
      Method = c("Biopsy alone", "PROMETheus model"),
      Sensitivity = c(biopsy_sensitivity, model_sensitivity),
      Specificity = c(biopsy_specificity, model_specificity),
      N_High_Risk = c(sum(actual_high), sum(actual_high)),
      stringsAsFactors = FALSE
    )
  )
  
  return(results)
}

# =============================================================================
# PART 6: execute
# =============================================================================

run_risk_group_analysis <- function() {
  
  cat("\n")
  cat("================================================================\n")
  cat("  RISK GROUP DIFFERENCES - EXTERNAL VALIDATION ANALYSIS\n")
  cat("================================================================\n\n")
  
  # load data
  loaded <- load_validation_data()
  validation_data <- loaded$data
  post <- loaded$post
  std_params <- loaded$std_params
  
  # compute risk classifications
  risk_data <- compute_risk_classifications(validation_data, post, std_params)
  
  # calculate metrics
  cat("\n--- Calculating Metrics ---\n")
  
  # overall accuracy
  overall_results <- calculate_overall_accuracy(risk_data)
  
  cat("\nOVERALL ACCURACY:\n")
  print(overall_results$overall)
  
  cat("\nHIGH-RISK DETECTION:\n")
  print(overall_results$high_risk)
  
  # stratified by risk group
  risk_stratified <- calculate_stratified_accuracy(risk_data)
  cat("\nACCURACY BY RISK GROUP:\n")
  print(risk_stratified)
  
  # stratified by surface
  surface_stratified <- calculate_surface_stratified_accuracy(risk_data)
  cat("\nACCURACY BY BIOPSY SURFACE:\n")
  print(surface_stratified)
  
  # plots
  cat("\n--- Creating Plots ---\n")
  plot_accuracy_by_risk_group(risk_stratified, save_plot = TRUE)
  plot_accuracy_by_surface(surface_stratified, save_plot = TRUE)
  plot_improvement_by_surface(surface_stratified, save_plot = TRUE)
  
  cat("\n================================================================\n")
  cat("  ANALYSIS COMPLETE\n")
  cat("================================================================\n")
  cat("\nOutput files saved to: validation_output1/figures/\n")
  
  # results for further use
  return(list(
    risk_data = risk_data,
    overall_results = overall_results,
    risk_stratified = risk_stratified,
    surface_stratified = surface_stratified
  ))
}

# =============================================================================
# run script
# =============================================================================

results <- run_risk_group_analysis()
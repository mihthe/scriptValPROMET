getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# risk function from app file 
risk_str <- function(size, mic, site) {
  
  size <- size/10 # in the simulation size is rnorm
  
  if (site == 4L) {
    
    if (mic <= 5) {
      if (size <= 2) {
        return("none")
      } else if (size <= 5) {
        return("very low")
      } else if (size <= 10) {
        return("low") 
      } else {
        return("moderate")
      }
      
    } else { # mic > 5
      if (size <= 2) {
        return("none")
      } else if (size <= 5) {
        return("moderate")
      } else if (size <= 10) {
        return("high") 
      } else {
        return("high")
      }
    }
    
  } else { 
    
    if (mic <= 5) {
      if (size <= 2) {
        return("none")
      } else if (size <= 5) {
        return("low")
      } else if (size <= 10) {
        return("moderate") 
      } else {
        return("high")
      }
      
    } else { # mic > 5
      if (size <= 2) {
        return("high")
      } else if (size <= 5) {
        return("high")
      } else if (size <= 10) {
        return("high") 
      } else {
        return("high")
      }
    }
  }
}

# Vectorized version for applying to entire dataframe or vector
risk_str_vectorized <- function(size, mic, site) {
  
  # handle single values and vectors
  n <- max(length(size), length(mic), length(site))
  
  # initialize output vector
  risk <- character(n)
  
  # risk function to each case
  for (i in 1:n) {
    risk[i] <- risk_str(
      size = size[i], 
      mic = mic[i], 
      site = site[i]
    )
  }
  
  return(risk)
}

# create risk plots file

# libraries
#library(rethinking)


###########################################
# FUNCTION 1: Compute Risk Classifications
###########################################

# adapting from app file
compute_risk_classifications <- function(db, post, use_predictions = TRUE) {
  # Computes risk classes based on biopsy and surgery mitotic counts
  #   check data frame has columns Si (size), m_bio, m_sur, Su, L
  #   posterior samples from the fitted model
  #   use_predictions: if TRUE, uses model predictions for surgery; if FALSE, uses observed m_sur
  #
  # Returns:
  #   data frame with original data plus risk classifications
  
  n_cases <- nrow(db) # using validation from simulation file
  
  # 1. risk based on BIOPSY mitotic count
  risk_biopsy <- risk_str_vectorized(
    size = db$Size, 
    mic = db$m_bio, 
    site = db$Site
  )
  
  # 2. risk based on SURGERY mitotic count 
  risk_surgery_actual <- risk_str_vectorized(
    size = db$Size, 
    mic = db$m_sur, 
    site = db$Site
  )
  
  # 3. risk based on PREDICTED surgery mitotic count
  if (use_predictions) {
    
    # link function to predict lambda (expected mitotic count)  
    link <- function(Si, m_bio, Su, L) {
      # std w\ training data parameters
      Si <- (Si - 62.175) / 48.07636
      m_bio <- (m_bio - 3.6625) / 8.066094
      Su <- (Su - 14.51) / 9.335991
      
      mu <- with(post, {
        a + b[, L] * Si + g * m_bio + e[, L] * Su
      })
      lambda <- exp(mu)
      return(lambda)
       #sim_m_sur <- rpois(4000, lambda)
       #return(sim_m_sur)
    }
    
    # predictions for each case
    risk_surgery_predicted <- character(n_cases)
    predicted_mitosis_m <- numeric(n_cases)
    
    for (i in 1:n_cases) {
      # lambda distribution for this case
      sim_m_sur <- link(
        Si = db$Size[i], 
        m_bio = db$m_bio[i], 
        Su = db$Su[i], 
        L = db$Site[i]
      )
      
      # median of lambda as predicted mitotic count
      pred_mitosis <- median(sim_m_sur)
      predicted_mitosis_m[i] <- pred_mitosis
      
      # risk class from predicted mitosis
      risk_surgery_predicted[i] <- risk_str(
        size = db$Size[i], 
        mic = round(pred_mitosis), 
        site = db$Site[i]
      )
      
      if (i %% 100 == 0) cat("  Processed", i, "cases\n")
    }
    
  } else {
    risk_surgery_predicted <- risk_surgery_actual
    predicted_mitosis_m <- db$m_sur
  }
  
  # everything into output data frame
  result <- data.frame(
    case_id = 1:n_cases,
    size = db$Site,
    location = db$Size,
    m_bio = db$m_bio,
    m_sur = db$m_sur,
    surface_bio = db$Su,
    risk_biopsy = risk_biopsy,
    risk_surgery_actual = risk_surgery_actual,
    risk_surgery_predicted = risk_surgery_predicted,
    predicted_mitosis = predicted_mitosis_m,
    stringsAsFactors = FALSE
  )
  
  return(result)
}


#########################################
# FUNCTION 2: Calculate Accuracy Metrics
#########################################

calculate_accuracy_metrics <- function(risk_data) {
  # comparing predictions to actual outcomes
  #   risk_data is output from compute_risk_classifications()
  #
  # Returns:
  #   list with accuracy metrics and confusion matrices
  
  # risk levels in order
  risk_levels <- c("none", "very low", "low", "moderate", "high")
  
  # convert to factors with ordered levels
  actual <- factor(risk_data$risk_surgery_actual, levels = risk_levels)
  predicted <- factor(risk_data$risk_surgery_predicted, levels = risk_levels)
  biopsy <- factor(risk_data$risk_biopsy, levels = risk_levels)
  
  # 1. Overall accuracy: Model prediction vs. Actual surgery
  correct_model <- sum(predicted == actual, na.rm = TRUE)
  accuracy_model <- correct_model / nrow(risk_data)
  
  # 2. Baseline accuracy: Biopsy vs. Actual surgery
  correct_biopsy <- sum(biopsy == actual, na.rm = TRUE)
  accuracy_biopsy <- correct_biopsy / nrow(risk_data)
  
  # 3. Confusion matrix: Model prediction vs. Actual
  conf_matrix_model <- table(
    Predicted = predicted, 
    Actual = actual, 
    useNA = "no"
  )
  
  # 4. Confusion matrix: Biopsy vs. Actual
  conf_matrix_biopsy <- table(
    Biopsy = biopsy, 
    Actual = actual, 
    useNA = "no"
  )
  
  # 5. sensitivity and specificity for "high risk" category
  high_risk_actual <- actual == "high"
  high_risk_predicted <- predicted == "high"
  high_risk_biopsy <- biopsy == "high"
  
  # Model performance for high risk
  sensitivity_model <- sum(high_risk_predicted & high_risk_actual, na.rm = TRUE) / 
    sum(high_risk_actual, na.rm = TRUE)
  specificity_model <- sum(!high_risk_predicted & !high_risk_actual, na.rm = TRUE) / 
    sum(!high_risk_actual, na.rm = TRUE)
  
  # Biopsy performance for high risk
  sensitivity_biopsy <- sum(high_risk_biopsy & high_risk_actual, na.rm = TRUE) / 
    sum(high_risk_actual, na.rm = TRUE)
  specificity_biopsy <- sum(!high_risk_biopsy & !high_risk_actual, na.rm = TRUE) / 
    sum(!high_risk_actual, na.rm = TRUE)
  
  # 6. Risk category shifts
  risk_increase <- sum(as.numeric(predicted) > as.numeric(biopsy), na.rm = TRUE)
  risk_decrease <- sum(as.numeric(predicted) < as.numeric(biopsy), na.rm = TRUE)
  risk_unchanged <- sum(as.numeric(predicted) == as.numeric(biopsy), na.rm = TRUE)
  
  # results
  results <- list(
    overall_accuracy = data.frame(
      Method = c("Biopsy alone", "Model prediction"),
      Accuracy = c(accuracy_biopsy, accuracy_model),
      N_correct = c(correct_biopsy, correct_model),
      N_total = c(nrow(risk_data), nrow(risk_data))
    ),
    
    high_risk_performance = data.frame(
      Method = c("Biopsy alone", "Model prediction"),
      Sensitivity = c(sensitivity_biopsy, sensitivity_model),
      Specificity = c(specificity_biopsy, specificity_model)
    ),
    
    risk_changes = data.frame(
      Change = c("Increased", "Decreased", "Unchanged"),
      N = c(risk_increase, risk_decrease, risk_unchanged),
      Proportion = c(risk_increase, risk_decrease, risk_unchanged) / nrow(risk_data)
    ),
    
    confusion_matrix_model = conf_matrix_model,
    confusion_matrix_biopsy = conf_matrix_biopsy,
    
    detailed_data = risk_data
  )
  
  return(results)
}


###########################################
# FUNCTION 3: Plot Accuracy Results 
###########################################

plot_accuracy_results <- function(accuracy_results, save_plots = TRUE, 
                                  output_dir = "output/figures/") {
  
  #   accuracy_results is output from calculate_accuracy_metrics()
  #   save_plots: if TRUE, saves plots as JPEG files
  #   "output/figures/": directory to save plots
  
  if (save_plots & !dir.exists("output/figures/")) {
    dir.create("output/figures/", recursive = TRUE)
  }
  
  risk_levels <- c("none", "very low", "low", "moderate", "high")
  risk_colors <- c("darkgreen", "green3", "gold", "orange", "red")
  
  # ============================================================
  # PLOT 1: Overall Accuracy Comparison
  # ============================================================
  
  if (save_plots) {
    jpeg(paste0("output/figures/", "accuracy_comparison_validation.jpg"), 
         units = "in", width = 8, height = 6, res = 300)
  }
  
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

  acc_data <- accuracy_results$overall_accuracy
  
  barplot(
    acc_data$Accuracy * 100,
    names.arg = acc_data$Method,
    ylim = c(0, 100),
    ylab = "Accuracy (%)",
    main = "Risk Classification Accuracy",
    col = c("skyblue", "seagreen3"),
    border = NA,
    las = 1
  )
  
  # percentage labels on bars
  text(
    x = c(0.7, 1.9),
    y = acc_data$Accuracy * 100 + 3,
    labels = paste0(round(acc_data$Accuracy * 100, 1), "%"),
    font = 2
  )
  
  # sample size
  mtext(paste("N =", acc_data$N_total[1]), side = 3, line = 0.5, cex = 0.9)
  
  if (save_plots) dev.off()
  
  
  # ============================================================
  # PLOT 2: Confusion Matrix Heatmap - Model Predictions
  # ============================================================
  
  if (save_plots) {
    jpeg(paste0("output/figures/", "confusion_matrix_model_validation.jpg"), 
         units = "in", width = 6, height = 6, res = 300)
  }
  
  par(mfrow = c(1, 1), 
      mar = c(9, 9, 5, 6),      
      mgp = c(6, 1, 0))
  
  conf_mat <- accuracy_results$confusion_matrix_model
  
  # check if confusion matrix has data
  if (sum(conf_mat) == 0) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
         main = "Confusion Matrix: No Data")
    text(1, 1, "No predictions available", cex = 1.5)
    if (save_plots) dev.off()
  } else {
    
    # normalize by column for proportions
    conf_mat_prop <- prop.table(conf_mat, margin = 2)
    conf_mat_prop[is.nan(conf_mat_prop)] <- 0
    
    # NaN with 0 (happens when a column has all zeros)
    conf_mat_prop[is.nan(conf_mat_prop)] <- 0
    
    # color palette
    n_colors <- 100
    color_palette <- colorRampPalette(c("white", "steelblue", "darkblue"))(n_colors)
    
    # max value
    max_prop <- max(conf_mat_prop, na.rm = TRUE)
    if (!is.finite(max_prop) || max_prop == 0) max_prop <- 1
    
    # heatmap
    image(
      x = 1:ncol(conf_mat_prop),
      y = 1:nrow(conf_mat_prop),
      z = t(conf_mat_prop[nrow(conf_mat_prop):1, ]),
      col = color_palette,
      xlab = "",                 # add with mtext
      ylab = "",                 # add with mtext
      main = "Confusion Matrix: Model Predictions",
      axes = FALSE,
      zlim = c(0, max_prop)
    )
    
    # axis titles manually with precise positioning
    mtext("Actual Risk (Surgery)", side = 1, line = 7, cex = 1.1)
    mtext("Predicted Risk (Model)", side = 2, line = 6, cex = 1.1)
    
    # axes
    axis(1, at = 1:ncol(conf_mat_prop), 
         labels = colnames(conf_mat_prop), 
         las = 2,                
         cex.axis = 0.9)         
    
    axis(2, at = 1:nrow(conf_mat_prop), 
         labels = rev(rownames(conf_mat_prop)), 
         las = 1,                
         cex.axis = 0.9)
    
    # grid
    abline(h = 0.5:(nrow(conf_mat_prop) + 0.5), col = "white", lwd = 2)
    abline(v = 0.5:(ncol(conf_mat_prop) + 0.5), col = "white", lwd = 2)
    
    # text with counts and proportions
    for (i in 1:nrow(conf_mat)) {
      for (j in 1:ncol(conf_mat)) {
        count <- conf_mat[i, j]
        prop <- conf_mat_prop[i, j]
        
        if (count > 0 && is.finite(prop)) {
          text(
            x = j,
            y = nrow(conf_mat) - i + 1,
            labels = paste0(count, "\n(", round(prop * 100, 1), "%)"),
            col = if(prop > 0.5) "white" else "black",
            font = 2,
            cex = 0.8
          )
        }
      }
    }
    
    # # color legend
    # legend_vals <- seq(0, max_prop, length.out = 5)
    # legend(
    #   x = ncol(conf_mat_prop) + 0.8,
    #   y = nrow(conf_mat_prop),
    #   legend = paste0(round(legend_vals * 100, 0), "%"),
    #   fill = color_palette[seq(1, n_colors, length.out = 5)],
    #   title = "Proportion",
    #   xpd = TRUE,
    #   bty = "n",
    #   cex = 0.9
    # )
    
    if (save_plots) dev.off()
  }
  
  
  # ============================================================
  # PLOT 3: Confusion Matrix Heatmap - Biopsy Only
  # ============================================================
  
  if (save_plots) {
    jpeg(paste0("output/figures/", "confusion_matrix_biopsy_validation.jpg"), 
         units = "in", width = 6, height = 6, res = 300)
  }
  
  par(mfrow = c(1, 1), 
      mar = c(9, 8, 5, 4),      
      mgp = c(5, 1, 0))  
  
  conf_mat_bio <- accuracy_results$confusion_matrix_biopsy
  
  # check if confusion matrix has data
  if (sum(conf_mat_bio) == 0) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
         main = "Confusion Matrix: No Data")
    text(1, 1, "No biopsy data available", cex = 1.5)
    if (save_plots) dev.off()
  } else {
    
    conf_mat_bio_prop <- prop.table(conf_mat_bio, margin = 2)
    conf_mat_bio_prop[is.nan(conf_mat_bio_prop)] <- 0
    
    # color palette
    n_colors <- 100
    color_palette <- colorRampPalette(c("white", "steelblue", "darkblue"))(n_colors)
    
    # max value 
    max_prop_bio <- max(conf_mat_bio_prop, na.rm = TRUE)
    if (!is.finite(max_prop_bio) || max_prop_bio == 0) max_prop_bio <- 1
    
    image(
      x = 1:ncol(conf_mat_bio_prop),
      y = 1:nrow(conf_mat_bio_prop),
      z = t(conf_mat_bio_prop[nrow(conf_mat_bio_prop):1, ]),
      col = color_palette,
      xlab = "Actual Risk (Surgery)",
      ylab = "Biopsy Risk",
      main = "Confusion Matrix: Biopsy Alone",
      axes = FALSE,
      zlim = c(0, max_prop_bio)
    )
    
    axis(1, at = 1:ncol(conf_mat_bio_prop), 
         labels = colnames(conf_mat_bio_prop), 
         las = 2,
         cex.axis = 0.9)           # Slightly smaller labels if needed
    
    axis(2, at = 1:nrow(conf_mat_bio_prop), 
         labels = rev(rownames(conf_mat_bio_prop)), 
         las = 1, 
         cex.axis = 0.9)
    
    abline(h = 0.5:(nrow(conf_mat_bio_prop) + 0.5), col = "white", lwd = 2)
    abline(v = 0.5:(ncol(conf_mat_bio_prop) + 0.5), col = "white", lwd = 2)
    
    for (i in 1:nrow(conf_mat_bio)) {
      for (j in 1:ncol(conf_mat_bio)) {
        count <- conf_mat_bio[i, j]
        prop <- conf_mat_bio_prop[i, j]
        
        if (count > 0 && is.finite(prop)) {
          text(
            x = j,
            y = nrow(conf_mat_bio) - i + 1,
            labels = paste0(count, "\n(", round(prop * 100, 1), "%)"),
            col = if(prop > 0.5) "white" else "black",
            font = 2,
            cex = 0.8
          )
        }
      }
    }
    
    legend_vals_bio <- seq(0, max_prop_bio, length.out = 5)
    legend(
      x = ncol(conf_mat_bio_prop) + 1.5,
      y = nrow(conf_mat_bio_prop),
      legend = paste0(round(legend_vals_bio * 100, 0), "%"),
      fill = color_palette[seq(1, n_colors, length.out = 5)],
      title = "Proportion",
      xpd = TRUE,
      bty = "n"
    )
    
    if (save_plots) dev.off()
  }
  
  
  # ============================================================
  # PLOT 4: Risk Category Changes (Biopsy -> Model Prediction)
  # ============================================================
  
  if (save_plots) {
    jpeg(paste0("output/figures/", "risk_category_changes_validation.jpg"), 
         units = "in", width = 8, height = 6, res = 300)
  }
  
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

  change_data <- accuracy_results$risk_changes
  
  barplot(
    change_data$Proportion * 100,
    names.arg = change_data$Change,
    ylim = c(0, max(change_data$Proportion * 100) * 1.2),
    ylab = "Percentage of Cases (%)",
    main = "Risk Category Changes: Biopsy to Model Prediction",
    col = c("coral", "steelblue", "gray70"),
    border = NA,
    las = 1
  )
  
  # labels
  text(
    x = c(0.7, 1.9, 3.1),
    y = change_data$Proportion * 100 + max(change_data$Proportion * 100) * 0.05,
    labels = paste0(
      round(change_data$Proportion * 100, 1), "%\n(n=", 
      change_data$N, ")"
    ),
    font = 2
  )
  
  if (save_plots) dev.off()
  
  
  # ============================================================
  # PLOT 5: Sensitivity & Specificity for High Risk
  # ============================================================
  
  if (save_plots) {
    jpeg(paste0("output/figures/", "high_risk_performance_validation.jpg"), 
         units = "in", width = 9, height = 7, res = 300)
  }
  
  par(mfrow = c(1, 1), 
      mar = c(6, 5, 5, 8),      
      xpd = FALSE)

  perf_data <- accuracy_results$high_risk_performance
  
  # matrix for grouped barplot
  perf_matrix <- t(as.matrix(perf_data[, c("Sensitivity", "Specificity")]))
  colnames(perf_matrix) <- perf_data$Method
  
  x_pos <- barplot(
    perf_matrix * 100,
    beside = TRUE,
    ylim = c(0, 115),            # extended to fit labels above 100%
    ylab = "Performance (%)",
    main = "Performance for 'High Risk' Classification",
    col = c("tomato", "skyblue"),
    border = NA,
    las = 1,                     # horizontal x-axis labels (option 1)
    # las = 2,                   # for rotated labels (option 2)
    names.arg = perf_data$Method,
    cex.names = 0.9
  )
  
  # percentage labels - positioned safely above bars
  for (i in 1:nrow(perf_matrix)) {
    for (j in 1:ncol(perf_matrix)) {
      label_y <- perf_matrix[i, j] * 100 + 5  # 5% above bar
      text(
        x = x_pos[i, j],
        y = label_y,
        labels = paste0(round(perf_matrix[i, j] * 100, 1), "%"),
        font = 2,
        cex = 0.9
      )
    }
  }
  
  # legend OUTSIDE the plot area
  par(xpd = TRUE)  # Allow drawing outside plot region
  legend(
    x = max(x_pos) + 1,          # Position to the right of bars
    
    y = 100,                      # Top of plot
    legend = c("Sensitivity", "Specificity"),
    fill = c("tomato", "skyblue"),
    border = NA,
    bty = "n",                   # No box around legend
    cex = 1
  )
  par(xpd = FALSE)  # Reset to default
  
  if (save_plots) dev.off()
  
  
  cat("\n========================================\n")
  cat("PLOTS CREATED SUCCESSFULLY\n")
  cat("========================================\n")
  if (save_plots) {
    cat("Saved to:", "output/figures/", "\n")
    cat("  - accuracy_comparison.jpg\n")
    cat("  - confusion_matrix_model.jpg\n")
    cat("  - confusion_matrix_biopsy.jpg\n")
    cat("  - risk_category_changes.jpg\n")
    cat("  - high_risk_performance.jpg\n")
  }
  
  invisible(NULL)
}

# get metrics

# computing risk classification
risk_data <- compute_risk_classifications(
  db = db,    # Pass the dataset
  post = post,            # Pass the posterior
  use_predictions = TRUE  # Use model predictions
)

# computing accuracy metrics
accuracy_results <- calculate_accuracy_metrics(risk_data)

# print overall accuracy results
print(accuracy_results$overall_accuracy)

# print for high risk perfomance
print(accuracy_results$high_risk_performance)

# print risk category changes
print(accuracy_results$risk_changes)

# print confusion matrix - model prediction
print(accuracy_results$confusion_matrix_model)

# print confusion matrix - biopsy
print(accuracy_results$confusion_matrix_biopsy)

# create and save plots
plot_accuracy_results(
  accuracy_results, 
  save_plots = TRUE,
  output_dir = "output/figures/"
)

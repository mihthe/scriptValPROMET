# see difference among groups
# consider cases with less or none biospy material

# library
library(rethinking)

source("risk_accuracy_assessment.R")

# run data simulation file

# run data analysis file 

# risk function from app file 
risk_str <- function(size, mic, site) {
  
  size <- size / 10
  
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

# vectorized version for applying to entire dataframe or vector
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

# adapting from app file to compute risk classification
compute_risk_classifications <- function(slim_sim, post, use_predictions = TRUE) {
  # risk classes based on biopsy and surgery mitotic counts
  # check data frame has columns Si (size), m_bio, m_sur, Su, L
  # posterior samples from the fitted model
  # use_predictions: if TRUE, uses model predictions for surgery; if FALSE, uses observed m_sur
  # output is data frame with original data plus risk classifications
  
  n_cases <- nrow(slim_sim) # using slim_sim from simulation file
  
  # 1. risk based on BIOPSY mitotic count
  risk_biopsy <- risk_str_vectorized(
    size = slim_sim$Si, 
    mic = slim_sim$m_bio, 
    site = slim_sim$L
  )
  
  # 2. risk based on SURGERY mitotic count 
  risk_surgery_actual <- risk_str_vectorized(
    size = slim_sim$Si, 
    mic = slim_sim$m_sur, 
    site = slim_sim$L
  )
  
  # 3. risk based on PREDICTED surgery mitotic count
  if (use_predictions) {
    
    # link function to predict lambda (expected mitotic count)
    m_link <- function(Si, m_bio, Su, L) {
      # std w\ training data parameters
      Si <- (Si - 62.175) / 48.07636
      m_bio <- (m_bio - 3.6625) / 8.066094
      Su <- (Su - 14.51) / 9.335991
      
      mu <- with(post, {
        a + b[, L] * Si + g * m_bio + e[, L] * Su
      })
      lambda <- exp(mu)
      return(lambda)
    }
    
    # predictions for each case
    risk_surgery_predicted <- character(n_cases)
    predicted_mitosis_median <- numeric(n_cases)
    
    for (i in 1:n_cases) {
      # lambda distribution for this case
      lambda <- m_link(
        Si = slim_sim$Si[i], 
        m_bio = slim_sim$m_bio[i], 
        Su = slim_sim$Su[i], 
        L = slim_sim$L[i]
      )
      
      # median of lambda as predicted mitotic count
      pred_mitosis <- median(lambda)
      predicted_mitosis_median[i] <- pred_mitosis
      
      # risk class from predicted mitosis
      risk_surgery_predicted[i] <- risk_str(
        size = slim_sim$Si[i], 
        mic = round(pred_mitosis), 
        site = slim_sim$L[i]
      )
      
      if (i %% 100 == 0) cat("  Processed", i, "cases\n")
    }
    
  } else {
    risk_surgery_predicted <- risk_surgery_actual
    predicted_mitosis_median <- slim_sim$m_sur
  }
  
  # everything into output data frame
  result <- data.frame(
    case_id = 1:n_cases,
    size = slim_sim$Si,
    location = slim_sim$L,
    m_bio = slim_sim$m_bio,
    m_sur = slim_sim$m_sur,
    surface_bio = slim_sim$Su,
    risk_biopsy = risk_biopsy,
    risk_surgery_actual = risk_surgery_actual,
    risk_surgery_predicted = risk_surgery_predicted,
    predicted_mitosis = predicted_mitosis_median,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# computing risk classification
risk_data <- compute_risk_classifications(
  slim_sim = slim_sim,    # Pass the dataset
  post = post,            # Pass the posterior
  use_predictions = TRUE  # Use model predictions
)

#######################################
# FUNCTION 1: Stratified Accuracy by Risk Group
#######################################

calculate_stratified_accuracy <- function(risk_data) {
  # for each actual risk category
  # risk_data is output from compute_risk_classifications()
  
  risk_levels <- c("none", "very low", "low", "moderate", "high")
  
  # initialize results
  results_list <- list()
  
  for (risk_level in risk_levels) {
    # cases with this actual risk
    subset_data <- risk_data[risk_data$risk_surgery_actual == risk_level, ]
    
    if (nrow(subset_data) > 0) {
      # model accuracy for this risk group
      correct_model <- sum(subset_data$risk_surgery_predicted == subset_data$risk_surgery_actual)
      accuracy_model <- correct_model / nrow(subset_data)
      
      # biopsy accuracy for this risk group
      correct_biopsy <- sum(subset_data$risk_biopsy == subset_data$risk_surgery_actual)
      accuracy_biopsy <- correct_biopsy / nrow(subset_data)
      
      # store results
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
  
  # all results
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  
  return(results)
}


#######################################
# FUNCTION 2: Stratified Accuracy by Biopsy Surface
#######################################

calculate_surface_stratified_accuracy <- function(risk_data) {
  # for different biopsy surface sizes
  # risk_data is output from compute_risk_classifications()
  
  # surface categories
  risk_data$surface_category <- cut(
    risk_data$surface_bio,
    breaks = c(0, 5, 10, 15, 20, 23.5),
    labels = c("0-5 HPF", "5-10 HPF", "10-15 HPF", "15-20 HPF", "20-23.5 HPF"),
    include.lowest = TRUE
  )
  
  # initialize results
  results_list <- list()
  
  surface_cats <- levels(risk_data$surface_category)
  
  for (cat in surface_cats) {
    # cases with this surface size
    subset_data <- risk_data[risk_data$surface_category == cat, ]
    
    if (nrow(subset_data) > 0) {
      # model accuracy for this surface category
      correct_model <- sum(subset_data$risk_surgery_predicted == subset_data$risk_surgery_actual, 
                           na.rm = TRUE)
      accuracy_model <- correct_model / nrow(subset_data)
      
      # biopsy accuracy for this surface category
      correct_biopsy <- sum(subset_data$risk_biopsy == subset_data$risk_surgery_actual, 
                            na.rm = TRUE)
      accuracy_biopsy <- correct_biopsy / nrow(subset_data)
      
      # mean predicted vs actual mitotic count
      mean_predicted <- mean(subset_data$predicted_mitosis, na.rm = TRUE)
      mean_actual <- mean(subset_data$m_sur, na.rm = TRUE)
      
      # store results
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
  
  # all results
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  
  return(results)
}


#######################################
# FUNCTION 3: Plot Accuracy by Risk Group
#######################################

plot_accuracy_by_risk_group <- function(stratified_results, save_plot = TRUE,
                                        output_dir = "output/figures/") {
  # plot accuracy by actual risk group
  
  if (save_plot) {
    jpeg(paste0(output_dir, "accuracy_by_risk_group.jpg"), 
         units = "in", width = 10, height = 6, res = 300)
  }
  
  par(mfrow = c(1, 1), mar = c(8, 5, 4, 2))
  
  # data for barplot
  acc_matrix <- t(as.matrix(stratified_results[, c("Biopsy_Accuracy", "Model_Accuracy")]))
  colnames(acc_matrix) <- stratified_results$Risk_Group
  
  # barplot
  bp <- barplot(
    acc_matrix * 100,
    beside = TRUE,
    ylim = c(0, 110),
    ylab = "Accuracy (%)",
    xlab = "",
    main = "Risk Classification Accuracy by Actual Risk Group",
    col = c("skyblue", "seagreen3"),
    border = NA,
    las = 2,
    legend.text = c("Biopsy alone", "Model prediction"),
    args.legend = list(x = "topright", bty = "n")
  )
  
  # percentage labels
  for (i in 1:nrow(acc_matrix)) {
    for (j in 1:ncol(acc_matrix)) {
      text(
        x = bp[i, j],
        y = acc_matrix[i, j] * 100 + 3,
        labels = paste0(round(acc_matrix[i, j] * 100, 1), "%"),
        font = 2,
        cex = 0.8
      )
    }
  }
  
  # sample sizes below bars
  mtext("Actual Risk Group (Surgery)", side = 1, line = 5, font = 2)
  for (j in 1:ncol(acc_matrix)) {
    text(
      x = mean(bp[, j]),
      y = -8,
      labels = paste0("n=", stratified_results$N_Cases[j]),
      cex = 0.8,
      xpd = TRUE
    )
  }
  
  if (save_plot) dev.off()
  
  invisible(NULL)
}


#######################################
# FUNCTION 4: Plot Accuracy by Biopsy Surface
#######################################

plot_accuracy_by_surface <- function(surface_results, save_plot = TRUE,
                                     output_dir = "output/figures/") {
  # plot accuracy by biopsy surface size
  
  if (save_plot) {
    jpeg(paste0(output_dir, "accuracy_by_surface.jpg"), 
         units = "in", width = 10, height = 6, res = 300)
  }
  
  par(mfrow = c(1, 1), mar = c(9, 5, 4, 2))
  
  # data for barplot
  acc_matrix <- t(as.matrix(surface_results[, c("Biopsy_Accuracy", "Model_Accuracy")]))
  colnames(acc_matrix) <- surface_results$Surface_Category
  
  # barplot
  bp <- barplot(
    acc_matrix * 100,
    beside = TRUE,
    ylim = c(0, 110),
    ylab = "Accuracy (%)",
    xlab = "",
    main = "Risk Classification Accuracy by Biopsy Surface Size",
    col = c("skyblue", "seagreen3"),
    border = NA,
    las = 2,
    legend.text = c("Biopsy alone", "Model prediction"),
    args.legend = list(x = "topright", bty = "n")
  )
  
  # percentage labels
  for (i in 1:nrow(acc_matrix)) {
    for (j in 1:ncol(acc_matrix)) {
      text(
        x = bp[i, j],
        y = acc_matrix[i, j] * 100 + 3,
        labels = paste0(round(acc_matrix[i, j] * 100, 1), "%"),
        font = 2,
        cex = 0.8
      )
    }
  }
  
  # sample sizes
  mtext("Biopsy Surface (High Power Fields)", side = 1, line = 6, font = 2)
  for (j in 1:ncol(acc_matrix)) {
    text(
      x = mean(bp[, j]),
      y = -8,
      labels = paste0("n=", surface_results$N_Cases[j]),
      cex = 0.8,
      xpd = TRUE
    )
  }
  
  if (save_plot) dev.off()
  
  invisible(NULL)
}


#######################################
# FUNCTION 5: Plot Improvement by Surface
#######################################

plot_improvement_by_surface <- function(surface_results, save_plot = TRUE,
                                        output_dir = "output/figures/") {
  # plot the improvement (model - biopsy accuracy) by surface size
  
  if (save_plot) {
    jpeg(paste0(output_dir, "improvement_by_surface.jpg"), 
         units = "in", width = 8, height = 6, res = 300)
  }
  
  par(mfrow = c(1, 1), mar = c(9, 5, 4, 2))
  
  # barplot of improvement
  bp <- barplot(
    surface_results$Improvement * 100,
    names.arg = surface_results$Surface_Category,
    ylim = c(-10, max(surface_results$Improvement * 100) + 10),
    ylab = "Accuracy Improvement (%)",
    xlab = "",
    main = "Model Improvement Over Biopsy Alone by Surface Size",
    col = ifelse(surface_results$Improvement > 0, "seagreen3", "coral"),
    border = NA,
    las = 2
  )
  
  # horizontal line at 0
  abline(h = 0, lty = 2, lwd = 2)
  
  # value labels
  text(
    x = bp,
    y = surface_results$Improvement * 100 + ifelse(surface_results$Improvement > 0, 2, -2),
    labels = paste0(ifelse(surface_results$Improvement > 0, "+", ""),
                    round(surface_results$Improvement * 100, 1), "%"),
    font = 2,
    cex = 0.9
  )
  
  # sample sizes
  mtext("Biopsy Surface (High Power Fields)", side = 1, line = 6, font = 2)
  for (j in 1:length(bp)) {
    text(
      x = bp[j],
      y = -8,
      labels = paste0("n=", surface_results$N_Cases[j]),
      cex = 0.8,
      xpd = TRUE
    )
  }
  
  if (save_plot) dev.off()
  
  invisible(NULL)
}


#########################################
# metrics and plots
#########################################

# compute risk classifications
risk_data <- compute_risk_classifications(
  slim_sim = slim_sim,
  post = post,
  use_predictions = TRUE
)

# overall accuracy metrics
accuracy_results <- calculate_accuracy_metrics(risk_data)

# accuracy by risk group
risk_stratified <- calculate_stratified_accuracy(risk_data)
print(risk_stratified)

# accuracy by biopsy surface
surface_stratified <- calculate_surface_stratified_accuracy(risk_data)
print(surface_stratified)

# overall accuracy
print(accuracy_results$overall_accuracy)

# high risk performance
print(accuracy_results$high_risk_performance)

# plots
plot_accuracy_results(
  accuracy_results, 
  save_plots = TRUE,
  output_dir = "output/figures/"
)
plot_accuracy_by_risk_group(risk_stratified, save_plot = TRUE)
plot_accuracy_by_surface(surface_stratified, save_plot = TRUE)
plot_improvement_by_surface(surface_stratified, save_plot = TRUE)

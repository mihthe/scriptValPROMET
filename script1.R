# generating a single simulated dataset called 
# N_simulations times by the main validation script.

# load libraries
library(rethinking)  # for inv_logit and statistical functions
library(MASS)        # for mvrnorm (multivariate normal)

# set seed for reproducibility (will be overridden in main validation script)
set.seed(2025)

#config parametres

CONFIG <- list(
  # sample size
  N_patients = 80,              # n patients in this dataset (changeable)
  
  # location settings  
  N_locations = 4,              # n tumor locations
  location_names = c("Colon-rectum", "Duodenum", "Small intestine", "Stomach"),
  location_probs = c(0.10, 0.10, 0.15, 0.65),  # Realistic GIST distribution
  
  # measurement parameters
  max_HPF = 23.5,               # max high-power fields in biopsy (5mm²)
  
  # prior parameters (these define "reasonable" parameter ranges)
  prior = list(
    alpha_mean = -1,            # log-scale intercept
    alpha_sd = 0.2,
    
    b_bar_mean = -1,            # avg size effect (log-scale)
    b_bar_sd = 0.2,
    
    e_bar_mean = -1,            # avg surface effect (log-scale)  
    e_bar_sd = 0.2,
    
    g_mean = -1,                # biopsy mitotic count effect (no therapy)
    g_sd = 0.2,
    
    d_mean = -1,                # biopsy mitotic count effect (with response)
    d_sd = 0.2,
    
    sigma_b_rate = 2,           # exp prior rate for between-location SD
    sigma_e_rate = 2
  ),
  
  # therapy simulation parameters
  therapy_cutoff_quantile = 0.5  # Top 50% of risk get neoadjuvant therapy
)


#------------------------------------------------------------------------------
# FUNCTION: SIMULATE ONE DATASET
#------------------------------------------------------------------------------

simulate_gist_dataset <- function(config = CONFIG, seed = NULL, verbose = TRUE) {
  #'
  #' simulate a complete GIST dataset with known ground truth
  #' 
  #' @param config List of configuration parameters
  #' @param seed Random seed (for reproducibility)
  #' @param verbose Print progress messages?
  #' 
  #' @return List containing:
  #'   - data: Data frame ready for model fitting
  #'   - truth: Named vector of true parameter values
  #'   - metadata: Additional information about the simulation
  
  if (!is.null(seed)) set.seed(seed)
  if (verbose) cat("Simulating dataset...\n")
  
  N <- config$N_patients
  H <- config$N_locations
  
  
  #----------------------------------------------------------------------------
  # STEP 1: draw true parameters from priors
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Drawing true parameters from priors...\n")
  
  # global parameters
  alpha <- rnorm(1, config$prior$alpha_mean, config$prior$alpha_sd)
  g <- rnorm(1, config$prior$g_mean, config$prior$g_sd)
  d <- rnorm(1, config$prior$d_mean, config$prior$d_sd)
  
  # hierarchical parameters - location-specific effects
  b_bar <- rnorm(1, config$prior$b_bar_mean, config$prior$b_bar_sd)
  e_bar <- rnorm(1, config$prior$e_bar_mean, config$prior$e_bar_sd)
  
  sigma_b <- rexp(1, config$prior$sigma_b_rate)
  sigma_e <- rexp(1, config$prior$sigma_e_rate)
  
  # draw location-specific deviations
  b <- rnorm(H, b_bar, sigma_b)  # size effects by location
  e <- rnorm(H, e_bar, sigma_e)  # surface effects by location
  
  
  #----------------------------------------------------------------------------
  # STEP 2: simulate latent biological aggressiveness (U)
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Simulating latent aggressiveness...\n")
  
  # U is unmeasured biological factors
  # higher U → more aggressive tumor
  U <- rnorm(N, mean = 0, sd = 1)
  
  
  #----------------------------------------------------------------------------
  # STEP 3: simulate tumor locations
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Assigning tumor locations...\n")
  
  L <- sample(
    x = 1:H, 
    size = N, 
    replace = TRUE, 
    prob = config$location_probs
  )
  
  # check if all locations represented (important for hierarchical model)
  if (length(unique(L)) < H) {
    warning(paste("Only", length(unique(L)), "of", H, "locations represented in sample"))
  }
  
  
  #----------------------------------------------------------------------------
  # STEP 4: simulate tumor sizes (standardized)
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Simulating tumor sizes...\n")
  
  # size depends on location and biological aggressiveness
  
  # location-specific intercepts and slopes for size
  alpha_size <- rnorm(H, mean = 0, sd = 1)
  beta_size <- rnorm(H, mean = 0, sd = 1)
  
  # generate sizes (already on standardized scale)
  mu_size <- alpha_size[L] + beta_size[L] * U
  Size_std <- rnorm(N, mean = mu_size, sd = 1)
  
  
  #----------------------------------------------------------------------------
  # STEP 5: simulate biopsy surface area (standardized)
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Simulating biopsy surface areas...\n")
  
  # surface depends on location and size
  alpha_surf <- rnorm(H, mean = 0, sd = 1)
  beta_surf <- rnorm(H, mean = 0, sd = 1)
  
  # generate on logit scale, then transform
  mu_surf_logit <- 2 + alpha_surf[L] + beta_surf[L] * Size_std
  surf_logit <- rnorm(N, mean = mu_surf_logit, sd = 2)
  Surface_raw <- inv_logit(surf_logit) * config$max_HPF
  
  # standardize for model input
  Surface_std <- (Surface_raw - mean(Surface_raw)) / sd(Surface_raw)
  
  
  #----------------------------------------------------------------------------
  # STEP 6: simulate mitotic count on surgical specimen
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Simulating surgical specimen mitotic counts...\n")
  
  # this is what we're trying to predict!
  # log-linear model for Poisson rate parameter
  log_lambda_surg <- alpha + b[L] * Size_std + 
    log(config$max_HPF) # Offset for measurement area
  
  lambda_surg <- exp(log_lambda_surg)
  Mitosis_surg <- rpois(N, lambda = lambda_surg)
  
  
  #----------------------------------------------------------------------------
  # STEP 7: simulate mitotic count on biopsy
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Simulating biopsy mitotic counts...\n")
  
  # biopsy count reflects surgical count but with:
  # - Smaller surface area (sampling variability)
  # - Tumor heterogeneity
  
  log_lambda_bio <- log(Surface_raw) +  # offset for smaller biopsy area
    alpha + 
    b[L] * Size_std
  
  lambda_bio <- exp(log_lambda_bio)
  Mitosis_bio <- rpois(N, lambda = lambda_bio)
  
  # standardize for model input
  Mitosis_bio_std <- (Mitosis_bio - mean(Mitosis_bio)) / sd(Mitosis_bio)
  
  
  #----------------------------------------------------------------------------
  # STEP 8: simulate neoadjuvant therapy and response
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Simulating therapy decisions and responses...\n")
  
  # high-risk patients get therapy
  risk_score <- Size_std + Mitosis_bio_std
  risk_cutoff <- quantile(risk_score, config$therapy_cutoff_quantile)
  receives_therapy <- (risk_score > risk_cutoff)
  
  # response to therapy inversely related to aggressiveness
  # more aggressive tumors (high U) respond less
  prob_response <- 1 - inv_logit(U - 1.5)
  
  therapy_response <- rep(0, N)
  therapy_response[receives_therapy] <- rbinom(
    sum(receives_therapy), 
    size = 1, 
    prob = prob_response[receives_therapy]
  )
  
  # if therapy response, surgical mitotic count dramatically reduced
  Mitosis_surg[therapy_response == 1] <- rpois(
    sum(therapy_response == 1), 
    lambda = 0.1  # very low mitotic activity after successful therapy
  )
  
  
  #----------------------------------------------------------------------------
  # STEP 9: create indicator variables for model
  #----------------------------------------------------------------------------
  
  # n = 1 if no therapy response (use g coefficient)
  # y = 1 if therapy response (use d coefficient)
  n_indicator <- as.integer(therapy_response == 0)
  y_indicator <- as.integer(therapy_response == 1)
  
  
  #----------------------------------------------------------------------------
  # STEP 10: package results
  #----------------------------------------------------------------------------
  
  if (verbose) cat("  Packaging results...\n")
  
  # data list for Stan/ulam (matching your model structure)
  sim_data <- list(
    # observed data (what goes into the model)
    m_surg = Mitosis_surg,
    m_bio = Mitosis_bio_std,  # Standardized
    Si = Size_std,            # Standardized
    Su = Surface_std,         # Standardized
    L = L,
    n = n_indicator,
    y = y_indicator,
    
    # Metadata
    N = N
  )
  
  # Store TRUE parameter values (for SBC rank calculation)
  # NOTE: We only validate primary parameters, not hyperparameters (b_bar, e_bar)
  # because they have different distributions (marginal vs conditional)
  true_params <- c(
    alpha = alpha,
    b_1 = b[1],
    b_2 = b[2],
    b_3 = b[3],
    b_4 = b[4],
    e_1 = e[1],
    e_2 = e[2],
    e_3 = e[3],
    e_4 = e[4],
    g = g,
    d = d,
    sigma_b = sigma_b,
    sigma_e = sigma_e
  )
  
  # additional metadata
  metadata <- list(
    N_patients = N,
    N_locations = H,
    location_counts = table(L),
    therapy_count = sum(receives_therapy),
    response_count = sum(therapy_response),
    mean_mitosis_surg = mean(Mitosis_surg),
    mean_mitosis_bio = mean(Mitosis_bio),
    seed = seed
  )
  
  if (verbose) {
    cat("  Done! Summary:\n")
    cat(sprintf("    Patients: %d\n", N))
    cat(sprintf("    Locations: %s\n", 
                paste(table(L), collapse = ", ")))
    cat(sprintf("    Therapy: %d (%.1f%%)\n", 
                sum(receives_therapy), 100*mean(receives_therapy)))
    cat(sprintf("    Response: %d (%.1f%%)\n", 
                sum(therapy_response), 100*mean(therapy_response)))
  }
  
  return(list(
    data = sim_data,
    truth = true_params,
    metadata = metadata
  ))
}


#------------------------------------------------------------------------------
# FUNCTION: EXTRACT STANDARDIZATION PARAMETERS
#------------------------------------------------------------------------------

get_standardization_params <- function(sim_data) {
  #'
  #' Extract the standardization parameters used in simulation
  #' Needed for posterior predictive checks
  #' 
  #' @param sim_data Output from simulate_gist_dataset()$data
  #' @return List with means and SDs
  
  list(
    m_bio_mean = attr(sim_data$m_bio, 'scaled:center'),
    m_bio_sd = attr(sim_data$m_bio, 'scaled:scale'),
    Su_mean = attr(sim_data$Su, 'scaled:center'),
    Su_sd = attr(sim_data$Su, 'scaled:scale')
  )
}


#------------------------------------------------------------------------------
# TESTING CODE 
#------------------------------------------------------------------------------

if (interactive() && !exists("SOURCED_FROM_MAIN")) {
  
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("TESTING SIMULATION ENGINE\n")
  cat(rep("=", 78), "\n\n", sep = "")

  
  # Test 1: single simulation with default settings
  cat("Test 1: Simulating with default settings...\n")
  test_sim <- simulate_gist_dataset(verbose = TRUE)
  
  cat("\nData structure:\n")
  str(test_sim$data, max.level = 1)
  
  cat("\nTrue parameters:\n")
  print(round(test_sim$truth, 3))
  
  cat("\nMetadata:\n")
  print(test_sim$metadata)
  
  
  # Test 2: multiple simulations with different seeds
  cat("\n\nTest 2: Running 5 simulations with different seeds...\n")
  test_results <- lapply(1:5, function(i) {
    simulate_gist_dataset(seed = 1000 + i, verbose = FALSE)
  })
  
  cat("Success! All simulations completed.\n")
  
  # Test 3: check parameter variation
  cat("\nTest 3: Parameter variation across simulations:\n")
  alpha_values <- sapply(test_results, function(x) x$truth["alpha"])
  cat(sprintf("  alpha range: [%.2f, %.2f]\n", 
              min(alpha_values), max(alpha_values)))
  
  g_values <- sapply(test_results, function(x) x$truth["g"])
  cat(sprintf("  g range: [%.2f, %.2f]\n", 
              min(g_values), max(g_values)))
  
  cat("\n✓ simulation engine tests passed!\n\n")
}

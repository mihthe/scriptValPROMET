#if needed (commented out)
#install.packages(required_packages)

library(rstan)
library(loo)
library(bayesplot)
library(tidyverse)
library(dagitty)
library(ggdag)
library(posterior)

# set options for Stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# 1. optimal datasplitting

calculate_optimal_split <- function(n_parameters) {
  # based on Joseph (2022): γ* = 1/(√p + 1)
  gamma_star <- 1 / (sqrt(n_parameters) + 1)
  train_ratio <- 1 - gamma_star
  
  cat(sprintf("Number of parameters: %d\n", n_parameters))
  cat(sprintf("Optimal test ratio (γ*): %.3f\n", gamma_star))
  cat(sprintf("Optimal train:test split: %.1f:%.1f\n", 
              train_ratio * 100, gamma_star * 100))
  
  return(list(gamma = gamma_star, train_ratio = train_ratio))
}

# for PROMETheus: 4 β[L] + 4 γ[L] + α + δ + ε + σ_β + σ_γ = 13 parameters
optimal_split <- calculate_optimal_split(13)

# 2. datasplitting f

split_data_optimal <- function(data, gamma, location_var = "location", 
                               seed = 2024) {
  set.seed(seed)
  
  # stratified sampling by location to maintain distribution
  test_indices <- data %>%
    group_by(!!sym(location_var)) %>%
    slice_sample(prop = gamma) %>%
    pull(row_id)
  
  train_data <- data %>% filter(!row_id %in% test_indices)
  test_data <- data %>% filter(row_id %in% test_indices)
  
  # check balance
  cat("\nLocation distribution:\n")
  cat("Training set:\n")
  print(table(train_data[[location_var]]))
  cat("\nTest set:\n")
  print(table(test_data[[location_var]]))
  
  return(list(train = train_data, test = test_data))
}

# 3. simulating GIST data

simulate_gist_data <- function(n = 160, seed = 42) {
  set.seed(seed)
  
  # location: 1=colon-rectum, 2=duodenum, 3=small intestine, 4=stomach
  locations <- c("colon_rectum", "duodenum", "small_intestine", "stomach")
  location <- sample(1:4, n, replace = TRUE, prob = c(0.1, 0.15, 0.25, 0.5))
  
  # tumor dimension (size in mm)
  dimension <- rnorm(n, mean = 50, sd = 30)
  dimension <- pmax(dimension, 10)  # minimum 10mm
  
  # biopsy surface area (in HPF, 23.5 HPF = 5mm²)
  surface <- runif(n, min = 0.5, max = 23.5)
  
  # response to therapy (0 = no therapy or no response, 1 = response)
  therapy_prob <- rbinom(n, 1, 0.3)  # 30% received therapy
  response <- ifelse(therapy_prob == 1, rbinom(n, 1, 0.6), 0)
  
  # true parameters (from paper's posterior estimates)
  alpha <- 1.75
  
  # location-specific effects (β for dimension)
  beta <- c(3.81, 0.12, 0.30, 0.40)  # colon-rectum, duodenum, SI, stomach
  
  # location-specific effects (γ for surface)
  gamma <- c(0.38, -0.09, -1.22, -0.17)
  
  # biopsy count effect
  delta <- 1.01
  epsilon <- -0.90  # reduced effect when response to therapy
  
  # generate true mitotic count on biopsy
  lambda_biopsy <- exp(alpha - 0.5 + 
                         beta[location] * scale(dimension)[,1] * 0.3 +
                         rnorm(n, 0, 0.3))
  mitosis_biopsy <- rpois(n, lambda_biopsy)
  
  # generate mitotic count on surgical specimen
  lambda_specimen <- exp(alpha + 
                           beta[location] * scale(dimension)[,1] +
                           gamma[location] * scale(surface)[,1] +
                           ifelse(response == 0, 
                                  delta * log(mitosis_biopsy + 1),
                                  epsilon * log(mitosis_biopsy + 1)))
  
  mitosis_specimen <- rpois(n, lambda_specimen)
  
  # create dataframe
  data <- tibble(
    row_id = 1:n,
    location = factor(location, levels = 1:4, labels = locations),
    location_num = location,
    dimension = dimension,
    surface = surface,
    mitosis_biopsy = mitosis_biopsy,
    response = response,
    mitosis_specimen = mitosis_specimen
  )
  
  return(data)
}

# 4. causal DAG visualisation

create_prometheus_dag <- function() {
  #DAG
  dag <- dagify(
    MS ~ D + S + MB,
    MB ~ D + S,
    D ~ L,
    S ~ L,
    exposure = "MB",
    outcome = "MS",
    labels = c(
      MS = "Mitosis\nSpecimen",
      MB = "Mitosis\nBiopsy",
      D = "Dimension",
      S = "Surface",
      L = "Location"
    )
  )
  
  # plot
  p <- ggdag_status(dag, use_labels = "label", text = FALSE) +
    theme_dag() +
    labs(title = "PROMETheus Causal DAG",
         subtitle = "Directed Acyclic Graph for Mitotic Count Prediction")
  
  return(list(dag = dag, plot = p))
}

# 5. data for stan


prepare_stan_data <- function(data) {
  # standardize continuous predictors
  data <- data %>%
    mutate(
      dimension_std = as.numeric(scale(dimension)),
      surface_std = as.numeric(scale(surface)),
      log_mb_plus1 = log(mitosis_biopsy + 1)
    )
  
  # Stan data list
  stan_data <- list(
    N = nrow(data),
    N_loc = length(unique(data$location_num)),
    location = data$location_num,
    dimension = data$dimension_std,
    surface = data$surface_std,
    mitosis_biopsy = data$log_mb_plus1,
    response = data$response,
    mitosis_specimen = data$mitosis_specimen
  )
  
  # Store scaling parameters for inverse transformation
  attr(stan_data, "dimension_mean") <- mean(data$dimension)
  attr(stan_data, "dimension_sd") <- sd(data$dimension)
  attr(stan_data, "surface_mean") <- mean(data$surface)
  attr(stan_data, "surface_sd") <- sd(data$surface)
  
  return(stan_data)
}

# 6. stan model code


stan_model_code <- "
data {
  int<lower=0> N;              // number of observations
  int<lower=1> N_loc;          // number of locations
  array[N] int<lower=1, upper=N_loc> location;  // location index
  vector[N] dimension;         // standardized tumor dimension
  vector[N] surface;           // standardized biopsy surface
  vector[N] mitosis_biopsy;    // log(mitosis_biopsy + 1)
  array[N] int<lower=0, upper=1> response;  // therapy response indicator
  array[N] int<lower=0> mitosis_specimen;   // outcome
}

parameters {
  real alpha;                  // intercept
  vector[N_loc] beta_raw;      // location-specific dimension effects (raw)
  vector[N_loc] gamma_raw;     // location-specific surface effects (raw)
  real delta;                  // biopsy count effect (no response)
  real epsilon;                // biopsy count effect (with response)
  
  real<lower=0> sigma_beta;    // SD for beta hierarchy
  real<lower=0> sigma_gamma;   // SD for gamma hierarchy
}

transformed parameters {
  vector[N_loc] beta;          // location-specific dimension effects
  vector[N_loc] gamma;         // location-specific surface effects
  vector[N] lambda;            // expected mitotic count
  
  // Non-centered parameterization for better sampling
  beta = sigma_beta * beta_raw;
  gamma = sigma_gamma * gamma_raw;
  
  // Linear predictor
  for (i in 1:N) {
    lambda[i] = alpha + 
                beta[location[i]] * dimension[i] +
                gamma[location[i]] * surface[i] +
                (response[i] == 0 ? delta : epsilon) * mitosis_biopsy[i];
  }
}

model {
  // Priors
  alpha ~ normal(1.75, 0.5);
  beta_raw ~ std_normal();
  gamma_raw ~ std_normal();
  delta ~ normal(1, 0.3);
  epsilon ~ normal(-0.9, 0.3);
  
  sigma_beta ~ exponential(1);
  sigma_gamma ~ exponential(1);
  
  // Likelihood
  mitosis_specimen ~ poisson_log(lambda);
}

generated quantities {
  vector[N] log_lik;           // log-likelihood for LOO-CV
  array[N] int y_rep;          // posterior predictive samples
  
  for (i in 1:N) {
    log_lik[i] = poisson_log_lpmf(mitosis_specimen[i] | lambda[i]);
    y_rep[i] = poisson_log_rng(lambda[i]);
  }
}
"


# 7. save result until now in result_part1

main_part1 <- function() {
  cat("=== PROMETheus Model Validation - Part 1 ===\n\n")
  
  # Calculate optimal split
  cat("Step 1: Calculating optimal data split...\n")
  optimal <- calculate_optimal_split(13)
  
  # Simulate data
  cat("\nStep 2: Simulating GIST data...\n")
  gist_data <- simulate_gist_data(n = 160)
  cat(sprintf("Generated %d observations\n", nrow(gist_data)))
  
  # Visualize DAG
  cat("\nStep 3: Creating causal DAG...\n")
  dag_result <- create_prometheus_dag()
  print(dag_result$plot)
  
  # Split data
  cat("\nStep 4: Splitting data optimally...\n")
  split_data <- split_data_optimal(gist_data, optimal$gamma)
  
  # Prepare data for Stan
  cat("\nStep 5: Preparing data for Stan...\n")
  train_stan <- prepare_stan_data(split_data$train)
  test_stan <- prepare_stan_data(split_data$test)
  
  cat("\nPart 1 Complete! Data ready for model fitting.\n")
  cat("Proceed to Part 2 for model fitting and validation.\n")
  
  return(list(
    data = gist_data,
    split = split_data,
    train_stan = train_stan,
    test_stan = test_stan,
    stan_code = stan_model_code,
    dag = dag_result
  ))
}

# cxecute Part 1
 result_part1 <- main_part1()


# Model Fitting, Validation, and Diagnostics

# 1. prior predictive simulation

prior_predictive_check <- function(stan_code, stan_data, n_sim = 100) {
  cat("Running prior predictive simulation...\n")
  
  # Stan code sampling from prior only
  prior_code <- gsub(
    "mitosis_specimen ~ poisson_log\\(lambda\\);",
    "// mitosis_specimen ~ poisson_log(lambda);  // commented out for prior pred",
    stan_code
  )
  
  # compile and sample
  prior_model <- stan_model(model_code = prior_code)
  
  prior_fit <- sampling(
    prior_model,
    data = stan_data,
    chains = 2,
    iter = 1000,
    warmup = 500,
    refresh = 0
  )
  
  # extract prior predictive samples
  y_rep_prior <- extract(prior_fit)$y_rep
  
  # plot
  p <- ppc_dens_overlay(
    y = rep(NA, stan_data$N),  # no observed data
    yrep = y_rep_prior[1:min(50, n_sim), ]
  ) +
    labs(
      title = "Prior Predictive Distribution",
      subtitle = "Expected mitotic counts before seeing data",
      x = "Mitotic Count",
      y = "Density"
    ) +
    xlim(0, 100)
  
  print(p)
  
  # summary statistics
  prior_summary <- tibble(
    mean_lambda = mean(y_rep_prior),
    sd_lambda = sd(y_rep_prior),
    q05 = quantile(y_rep_prior, 0.05),
    q50 = quantile(y_rep_prior, 0.50),
    q95 = quantile(y_rep_prior, 0.95),
    prop_over_50 = mean(y_rep_prior > 50)
  )
  
  cat("\nPrior predictive summary:\n")
  print(prior_summary)
  
  return(list(samples = y_rep_prior, plot = p, summary = prior_summary))
}

# 2. fit model

fit_prometheus_model <- function(stan_code, stan_data, 
                                 chains = 4, iter = 4000, warmup = 2000) {
  cat("Compiling Stan model...\n")
  model <- stan_model(model_code = stan_code)
  
  cat("Fitting model...\n")
  fit <- sampling(
    model,
    data = stan_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    thin = 2,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 500
  )
  
  return(fit)
}

# 3. convergence

check_convergence <- function(fit) {
  cat("\n=== Convergence Diagnostics ===\n\n")
  
  # check Rhat
  rhat <- rhat(fit)
  cat("Rhat statistics:\n")
  cat(sprintf("  Max Rhat: %.4f\n", max(rhat, na.rm = TRUE)))
  cat(sprintf("  Parameters with Rhat > 1.01: %d\n", 
              sum(rhat > 1.01, na.rm = TRUE)))
  
  # check effective sample size
  ess_bulk <- ess_bulk(fit)
  ess_tail <- ess_tail(fit)
  
  cat("\nEffective Sample Size:\n")
  cat(sprintf("  Min ESS bulk: %.0f\n", min(ess_bulk, na.rm = TRUE)))
  cat(sprintf("  Min ESS tail: %.0f\n", min(ess_tail, na.rm = TRUE)))
  
  # check divergences
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  
  cat(sprintf("\nDivergent transitions: %d\n", divergences))
  
  # trace plots
  params_to_plot <- c("alpha", "delta", "epsilon", "sigma_beta", "sigma_gamma")
  p_trace <- mcmc_trace(fit, pars = params_to_plot) +
    labs(title = "Trace Plots for Key Parameters")
  
  # rank plots (trankplots)
  p_rank <- mcmc_rank_overlay(fit, pars = params_to_plot) +
    labs(title = "Rank Plots for Chain Mixing")
  
  print(p_trace)
  print(p_rank)
  
  # diagnostic summary
  diagnostics <- list(
    rhat = rhat,
    ess_bulk = ess_bulk,
    ess_tail = ess_tail,
    divergences = divergences,
    converged = max(rhat, na.rm = TRUE) < 1.01 && divergences == 0
  )
  
  return(diagnostics)
}

# 4. posterior predictive checks

posterior_predictive_check <- function(fit, stan_data) {
  cat("\n=== Posterior Predictive Checks ===\n\n")
  
  # extract posterior predictive samples
  y_rep <- extract(fit)$y_rep
  y <- stan_data$mitosis_specimen
  
  # 1. density overlay
  p1 <- ppc_dens_overlay(y, y_rep[1:100, ]) +
    labs(title = "Posterior Predictive Density",
         subtitle = "Observed vs. Predicted Mitotic Counts")
  
  # 2. distribution of statistics
  p2 <- ppc_stat(y, y_rep, stat = "mean") +
    labs(title = "Posterior Predictive Check: Mean")
  
  p3 <- ppc_stat(y, y_rep, stat = "sd") +
    labs(title = "Posterior Predictive Check: SD")
  
  # 3. empirical CDF
  p4 <- ppc_ecdf_overlay(y, y_rep[1:50, ]) +
    labs(title = "Empirical CDF Comparison")
  
  # 4. intervals
  p5 <- ppc_intervals(y, y_rep, x = seq_along(y), prob = 0.5, prob_outer = 0.9) +
    labs(title = "50% and 90% Prediction Intervals",
         x = "Observation Index", y = "Mitotic Count")
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  
  # Bayesian p-values
  T_obs <- mean(y)
  T_rep <- apply(y_rep, 1, mean)
  p_value_mean <- mean(T_rep >= T_obs)
  
  cat(sprintf("\nBayesian p-value (mean): %.3f\n", p_value_mean))
  cat("(Values close to 0 or 1 indicate poor fit)\n")
  
  return(list(
    plots = list(p1, p2, p3, p4, p5),
    p_value = p_value_mean
  ))
}

# 5. model comparison, WAIC and LOO

compute_information_criteria <- function(fit) {
  cat("\n=== Information Criteria ===\n\n")
  
  # extract log-likelihood
  log_lik <- extract_log_lik(fit, parameter_name = "log_lik")
  
  # compute WAIC
  waic_result <- waic(log_lik)
  cat("WAIC (Widely Applicable Information Criterion):\n")
  print(waic_result)
  
  # compute LOO-CV
  loo_result <- loo(log_lik)
  cat("\n\nLOO-CV (Leave-One-Out Cross-Validation):\n")
  print(loo_result)
  
  # check Pareto k diagnostic
  cat("\n\nPareto k diagnostic:\n")
  cat(sprintf("  Good (k < 0.5): %d\n", sum(loo_result$diagnostics$pareto_k < 0.5)))
  cat(sprintf("  Ok (0.5 ≤ k < 0.7): %d\n", 
              sum(loo_result$diagnostics$pareto_k >= 0.5 & 
                    loo_result$diagnostics$pareto_k < 0.7)))
  cat(sprintf("  Bad (k ≥ 0.7): %d\n", sum(loo_result$diagnostics$pareto_k >= 0.7)))
  
  # plot Pareto k
  p <- plot(loo_result, label_points = TRUE) +
    labs(title = "PSIS Diagnostic: Pareto k values")
  print(p)
  
  return(list(waic = waic_result, loo = loo_result))
}

# 6. OOS validation

validate_out_of_sample <- function(fit, train_data, test_data) {
  cat("\n=== Out-of-Sample Validation ===\n\n")
  
  # extract posterior samples
  posterior <- extract(fit)
  n_samples <- length(posterior$alpha)
  
  # make predictions for test set
  n_test <- test_data$N
  y_pred <- matrix(0, n_samples, n_test)
  
  for (i in 1:n_samples) {
    lambda_test <- posterior$alpha[i] +
      posterior$beta[i, test_data$location] * test_data$dimension +
      posterior$gamma[i, test_data$location] * test_data$surface +
      ifelse(test_data$response == 0,
             posterior$delta[i],
             posterior$epsilon[i]) * test_data$mitosis_biopsy
    
    y_pred[i, ] <- rpois(n_test, exp(lambda_test))
  }
  
  # compute predictive metrics
  y_test <- test_data$mitosis_specimen
  y_pred_mean <- colMeans(y_pred)
  y_pred_median <- apply(y_pred, 2, median)
  
  # mean Absolute Error
  mae <- mean(abs(y_test - y_pred_mean))
  
  # root Mean Squared Error
  rmse <- sqrt(mean((y_test - y_pred_mean)^2))
  
  # coverage of 90% prediction intervals
  y_pred_lower <- apply(y_pred, 2, quantile, probs = 0.05)
  y_pred_upper <- apply(y_pred, 2, quantile, probs = 0.95)
  coverage_90 <- mean(y_test >= y_pred_lower & y_test <= y_pred_upper)
  
  # expected log predictive density
  elpd <- mean(log(colMeans(dpois(rep(y_test, each = n_samples), 
                                  exp(y_pred)))))
  
  cat(sprintf("Mean Absolute Error: %.3f\n", mae))
  cat(sprintf("Root Mean Squared Error: %.3f\n", rmse))
  cat(sprintf("90%% Prediction Interval Coverage: %.3f\n", coverage_90))
  cat(sprintf("Expected Log Predictive Density: %.3f\n", elpd))
  
  # plot predictions vs observed
  pred_df <- tibble(
    observed = y_test,
    predicted = y_pred_mean,
    lower = y_pred_lower,
    upper = y_pred_upper,
    location = factor(test_data$location)
  )
  
  p1 <- ggplot(pred_df, aes(x = observed, y = predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_point(aes(color = location), alpha = 0.6) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = location), 
                  alpha = 0.3, width = 0.5) +
    labs(title = "Out-of-Sample Predictions",
         subtitle = "With 90% prediction intervals",
         x = "Observed Mitotic Count",
         y = "Predicted Mitotic Count") +
    theme_minimal()
  
  print(p1)
  
  # residual plot
  pred_df$residual <- pred_df$observed - pred_df$predicted
  
  p2 <- ggplot(pred_df, aes(x = predicted, y = residual)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(aes(color = location), alpha = 0.6) +
    labs(title = "Residual Plot",
         x = "Predicted Mitotic Count",
         y = "Residual (Observed - Predicted)") +
    theme_minimal()
  
  print(p2)
  
  metrics <- list(
    mae = mae,
    rmse = rmse,
    coverage_90 = coverage_90,
    elpd = elpd,
    predictions = pred_df
  )
  
  return(metrics)
}

# 7. parameters int

summarize_parameters <- function(fit) {
  cat("\n=== Parameter Estimates ===\n\n")
  
  # extract posterior samples
  posterior <- extract(fit)
  
  # main parameters
  params <- c("alpha", "delta", "epsilon", "sigma_beta", "sigma_gamma")
  
  summary_stats <- data.frame(
    parameter = params,
    mean = sapply(params, function(p) mean(posterior[[p]])),
    sd = sapply(params, function(p) sd(posterior[[p]])),
    q2.5 = sapply(params, function(p) quantile(posterior[[p]], 0.025)),
    q50 = sapply(params, function(p) quantile(posterior[[p]], 0.50)),
    q97.5 = sapply(params, function(p) quantile(posterior[[p]], 0.975))
  )
  
  print(summary_stats)
  
  # location-specific effects
  cat("\n\nLocation-specific dimension effects (β):\n")
  beta_summary <- data.frame(
    location = 1:ncol(posterior$beta),
    mean = colMeans(posterior$beta),
    sd = apply(posterior$beta, 2, sd),
    q2.5 = apply(posterior$beta, 2, quantile, 0.025),
    q97.5 = apply(posterior$beta, 2, quantile, 0.975)
  )
  print(beta_summary)
  
  cat("\n\nLocation-specific surface effects (γ):\n")
  gamma_summary <- data.frame(
    location = 1:ncol(posterior$gamma),
    mean = colMeans(posterior$gamma),
    sd = apply(posterior$gamma, 2, sd),
    q2.5 = apply(posterior$gamma, 2, quantile, 0.025),
    q97.5 = apply(posterior$gamma, 2, quantile, 0.975)
  )
  print(gamma_summary)
  
  # visualization
  p <- mcmc_areas(fit, pars = params, prob = 0.89, prob_outer = 0.95) +
    labs(title = "Posterior Distributions of Key Parameters",
         subtitle = "89% and 95% credible intervals")
  print(p)
  
  return(list(
    main = summary_stats,
    beta = beta_summary,
    gamma = gamma_summary
  ))
}


main_part2 <- function(result_part1) {
  cat("\n=== PROMETheus Model Validation - Part 2 ===\n\n")
  
  # prior predictive check
  cat("Step 1: Prior Predictive Simulation...\n")
  prior_check <- prior_predictive_check(
    result_part1$stan_code,
    result_part1$train_stan
  )
  
  # fit the model
  cat("\nStep 2: Fitting Model to Training Data...\n")
  fit <- fit_prometheus_model(
    result_part1$stan_code,
    result_part1$train_stan,
    chains = 4,
    iter = 4000,
    warmup = 2000
  )
  
  # check convergence
  cat("\nStep 3: Convergence Diagnostics...\n")
  diagnostics <- check_convergence(fit)
  
  if (!diagnostics$converged) {
    warning("Model did not converge properly! Consider increasing adapt_delta or iterations.")
  }
  
  # posterior predictive checks
  cat("\nStep 4: Posterior Predictive Checks...\n")
  ppc <- posterior_predictive_check(fit, result_part1$train_stan)
  
  # information criteria
  cat("\nStep 5: Computing Information Criteria...\n")
  ic <- compute_information_criteria(fit)
  
  # out-of-sample validation
  cat("\nStep 6: Out-of-Sample Validation...\n")
  oos_metrics <- validate_out_of_sample(
    fit,
    result_part1$train_stan,
    result_part1$test_stan
  )
  
  # parameter interpretation
  cat("\nStep 7: Parameter Estimation and Interpretation...\n")
  params <- summarize_parameters(fit)
  
  cat("\n\n=== Validation Complete ===\n")
  cat("\nModel Performance Summary:\n")
  cat(sprintf("  - Convergence: %s\n", 
              ifelse(diagnostics$converged, "PASSED", "FAILED")))
  cat(sprintf("  - Test Set MAE: %.3f\n", oos_metrics$mae))
  cat(sprintf("  - Test Set RMSE: %.3f\n", oos_metrics$rmse))
  cat(sprintf("  - 90%% Coverage: %.3f\n", oos_metrics$coverage_90))
  cat(sprintf("  - Expected log pred density: %.3f\n", oos_metrics$elpd))
  
  return(list(
    fit = fit,
    diagnostics = diagnostics,
    ppc = ppc,
    information_criteria = ic,
    oos_validation = oos_metrics,
    parameters = params
  ))
}

# save in result_part2
 result_part2 <- main_part2(result_part1)

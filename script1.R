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
  # ensure numeric location column
  if(!"location_num" %in% names(data)) {
    stop("Data must have 'location_num' column with integer location codes")
  }
  
  # standardize continuous predictors
  data <- data %>%
    mutate(
      dimension_std = as.numeric(scale(dimension)),
      surface_std = as.numeric(scale(surface)),
      log_mb_plus1 = log(mitosis_biopsy + 1)
    )
  
  # Stan data list - use location_num 
  stan_data <- list(
    N = nrow(data),
    N_loc = 4,  # there are 4 locations
    location = as.integer(data$location_num),  # to integer
    dimension = data$dimension_std,
    surface = data$surface_std,
    mitosis_biopsy = data$log_mb_plus1,
    response = as.integer(data$response),  # ensure this is integer
    mitosis_specimen = as.integer(data$mitosis_specimen)  # and this
  )
  
  # store scaling parameters
  attr(stan_data, "dimension_mean") <- mean(data$dimension)
  attr(stan_data, "dimension_sd") <- sd(data$dimension)
  attr(stan_data, "surface_mean") <- mean(data$surface)
  attr(stan_data, "surface_sd") <- sd(data$surface)
  
  return(stan_data)
}

# 6. Stan model code


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
  
  # calculate optimal split
  cat("Step 1: Calculating optimal data split...\n")
  optimal <- calculate_optimal_split(13)
  
  # simulate data
  cat("\nStep 2: Simulating GIST data...\n")
  gist_data <- simulate_gist_data(n = 160)
  cat(sprintf("Generated %d observations\n", nrow(gist_data)))
  
  # visualize DAG
  cat("\nStep 3: Creating causal DAG...\n")
  dag_result <- create_prometheus_dag()
  print(dag_result$plot)
  
  # split data
  cat("\nStep 4: Splitting data optimally...\n")
  split_data <- split_data_optimal(gist_data, optimal$gamma)
  
  # prepare data for Stan
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

# execute Part 1
result_part1 <- main_part1()


# Model Fitting, Validation, and Diagnostics

# 1. prior predictive simulation

prior_predictive_check <- function(stan_code, stan_data, n_sim = 100) {
  cat("Running prior predictive simulation...\n")
  
  # modified Stan code sampling from prior only
  # remove log_lik from generated quantities, not using likelihood
  prior_code <- "
data {
  int<lower=0> N;
  int<lower=1> N_loc;
  array[N] int<lower=1, upper=N_loc> location;
  vector[N] dimension;
  vector[N] surface;
  vector[N] mitosis_biopsy;
  array[N] int<lower=0, upper=1> response;
}

parameters {
  real alpha;
  vector[N_loc] beta_raw;
  vector[N_loc] gamma_raw;
  real delta;
  real epsilon;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_gamma;
}

transformed parameters {
  vector[N_loc] beta = sigma_beta * beta_raw;
  vector[N_loc] gamma = sigma_gamma * gamma_raw;
  vector[N] lambda;
  
  for (i in 1:N) {
    lambda[i] = alpha + 
                beta[location[i]] * dimension[i] +
                gamma[location[i]] * surface[i] +
                (response[i] == 0 ? delta : epsilon) * mitosis_biopsy[i];
  }
}

model {
  // Priors only - NO LIKELIHOOD
  alpha ~ normal(1.75, 0.5);
  beta_raw ~ std_normal();
  gamma_raw ~ std_normal();
  delta ~ normal(1, 0.3);
  epsilon ~ normal(-0.9, 0.3);
  sigma_beta ~ exponential(1);
  sigma_gamma ~ exponential(1);
  
  // Likelihood is commented out for prior predictive sampling
  // mitosis_specimen ~ poisson_log(lambda);
}

generated quantities {
  // Only generate prior predictive samples, no log_lik
  array[N] int y_rep;
  
  for (i in 1:N) {
    y_rep[i] = poisson_log_rng(lambda[i]);
  }
}
"

# compile and sample from prior
cat("Compiling prior predictive model...\n")
prior_model <- stan_model(model_code = prior_code)

cat("Sampling from prior...\n")
prior_fit <- sampling(
  prior_model,
  data = stan_data,
  chains = 2,
  iter = 1000,
  warmup = 500,
  refresh = 0,
  show_messages = FALSE
)

# extract prior predictive samples
y_rep_prior <- extract(prior_fit)$y_rep

# create visualization
p <- ppc_dens_overlay(
  y = stan_data$mitosis_specimen,
  yrep = y_rep_prior[1:min(50, nrow(y_rep_prior)), ]
) +
  labs(
    title = "Prior Predictive Distribution",
    subtitle = "Simulated data from priors (light) vs. Observed data (dark)",
    x = "Mitotic Count on Specimen",
    y = "Density"
  ) +
  xlim(0, max(100, max(stan_data$mitosis_specimen) + 10))

print(p)

# summary statistics
prior_summary <- tibble(
  mean = mean(y_rep_prior),
  sd = sd(y_rep_prior),
  min = min(y_rep_prior),
  q05 = quantile(y_rep_prior, 0.05),
  q25 = quantile(y_rep_prior, 0.25),
  q50 = quantile(y_rep_prior, 0.50),
  q75 = quantile(y_rep_prior, 0.75),
  q95 = quantile(y_rep_prior, 0.95),
  max = max(y_rep_prior),
  prop_over_50 = mean(y_rep_prior > 50),
  prop_zero = mean(y_rep_prior == 0)
)

cat("\nPrior predictive summary:\n")
print(prior_summary)

# check if priors are reasonable
cat("\nPrior reasonableness check:\n")
if (prior_summary$prop_over_50 > 0.1) {
  cat("  ⚠ Warning: >10% of prior samples have mitotic count > 50\n")
  cat("     (This might be too high - consider tighter priors)\n")
} else {
  cat("  ✓ Priors appear reasonable\n")
}

if (prior_summary$prop_zero > 0.5) {
  cat("  ⚠ Warning: >50% of prior samples are zero\n")
  cat("     (This might be too conservative)\n")
} else {
  cat("  ✓ Prior allows reasonable variation\n")
}

return(list(
  samples = y_rep_prior, 
  plot = p, 
  summary = prior_summary,
  fit = prior_fit
))
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
  
  # summary to extract Rhat and ESS
  fit_summary <- summary(fit)$summary
  
  # extract Rhat values
  rhat_vals <- fit_summary[, "Rhat"]
  
  cat("Rhat statistics:\n")
  cat(sprintf("  Max Rhat: %.4f\n", max(rhat_vals, na.rm = TRUE)))
  cat(sprintf("  Parameters with Rhat > 1.01: %d\n", 
              sum(rhat_vals > 1.01, na.rm = TRUE)))
  
  # extract effective sample sizes
  ess_bulk_vals <- fit_summary[, "n_eff"]
  
  cat("\nEffective Sample Size:\n")
  cat(sprintf("  Min ESS: %.0f\n", min(ess_bulk_vals, na.rm = TRUE)))
  
  # check divergences
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  
  cat(sprintf("\nDivergent transitions: %d\n", divergences))
  
  # check max treedepth
  max_treedepth_exceeded <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 12)))
  cat(sprintf("Max treedepth exceeded: %d\n", max_treedepth_exceeded))
  
  # trace plots for key parameters
  params_to_plot <- c("alpha", "delta", "epsilon", "sigma_beta", "sigma_gamma")
  
  cat("\nGenerating diagnostic plots...\n")
  
  p_trace <- mcmc_trace(fit, pars = params_to_plot) +
    labs(title = "Trace Plots for Key Parameters")
  
  print(p_trace)
  
  # rank plots
  p_rank <- mcmc_rank_overlay(fit, pars = params_to_plot) +
    labs(title = "Rank Plots for Chain Mixing")
  
  print(p_rank)
  
  # diagnostic summary
  diagnostics <- list(
    rhat = rhat_vals,
    ess = ess_bulk_vals,
    divergences = divergences,
    max_treedepth_exceeded = max_treedepth_exceeded,
    converged = max(rhat_vals, na.rm = TRUE) < 1.01 && divergences == 0
  )
  
  # print warnings
  if (divergences > 0) {
    cat("\n⚠ WARNING: Model has divergent transitions!\n")
    cat("   This suggests the posterior is difficult to explore.\n")
    cat("   The model may still be usable if divergences are < 1% of samples.\n")
  }
  
  if (max_treedepth_exceeded > 100) {
    cat("\n⚠ WARNING: Many transitions exceeded max treedepth!\n")
    cat("   Consider increasing max_treedepth or reparameterizing.\n")
  }
  
  return(diagnostics)
}

# 4. posterior predictive checks

posterior_predictive_check <- function(fit, stan_data) {
  cat("\n=== Posterior Predictive Checks ===\n\n")
  
  # extract posterior predictive samples - use rstan::extract explicitly
  posterior <- rstan::extract(fit)
  y_rep <- posterior$y_rep
  y <- stan_data$mitosis_specimen
  
  # 1. density overlay
  p1 <- ppc_dens_overlay(y, y_rep[1:100, ]) +
    labs(title = "Posterior Predictive Density",
         subtitle = "Observed vs. Predicted Mitotic Counts")
  
  print(p1)
  
  # 2. distribution of statistics
  p2 <- ppc_stat(y, y_rep, stat = "mean") +
    labs(title = "Posterior Predictive Check: Mean")
  
  print(p2)
  
  p3 <- ppc_stat(y, y_rep, stat = "sd") +
    labs(title = "Posterior Predictive Check: SD")
  
  print(p3)
  
  # 3. empirical CDF
  p4 <- ppc_ecdf_overlay(y, y_rep[1:50, ]) +
    labs(title = "Empirical CDF Comparison")
  
  print(p4)
  
  # 4. intervals
  p5 <- ppc_intervals(y, y_rep, x = seq_along(y), prob = 0.5, prob_outer = 0.9) +
    labs(title = "50% and 90% Prediction Intervals",
         x = "Observation Index", y = "Mitotic Count")
  
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
  
  # extract posterior samples - use rstan::extract explicitly
  posterior <- rstan::extract(fit)
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
  
  # mean absolute error
  mae <- mean(abs(y_test - y_pred_mean))
  
  # root mean squared error
  rmse <- sqrt(mean((y_test - y_pred_mean)^2))
  
  # coverage of 90% prediction intervals
  y_pred_lower <- apply(y_pred, 2, quantile, probs = 0.05)
  y_pred_upper <- apply(y_pred, 2, quantile, probs = 0.95)
  coverage_90 <- mean(y_test >= y_pred_lower & y_test <= y_pred_upper)
  
  # expected log predictive density (c)
  lpd <- numeric(n_test)
  for (j in 1:n_test) {
    # for each test point, average the probability across posterior samples
    lpd[j] <- log(mean(dpois(y_test[j], y_pred[, j])))
  }
  elpd <- sum(lpd)
  
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
    geom_point(aes(color = location), alpha = 0.6, size = 3) +
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
    geom_point(aes(color = location), alpha = 0.6, size = 3) +
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
  
  # extract posterior samples - use rstan::extract explicitly
  posterior <- rstan::extract(fit)
  
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
  
  # skip prior predictive check (?)
  cat("Step 1: Prior Predictive Simulation...\n")
  cat("(Skipped - proceeding directly to model fitting)\n")
  cat("Note: Priors are based on paper's posterior estimates:\n")
  cat("  - alpha ~ normal(1.75, 0.5)\n")
  cat("  - delta ~ normal(1, 0.3)\n")
  cat("  - epsilon ~ normal(-0.9, 0.3)\n\n")
  prior_check <- NULL
  
  # fit the model
  cat("Step 2: Fitting Model to Training Data...\n")
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
    prior_check = prior_check,  # Will be NULL
    diagnostics = diagnostics,
    ppc = ppc,
    information_criteria = ic,
    oos_validation = oos_metrics,
    parameters = params
  ))
}

# save in result_part2
result_part2 <- main_part2(result_part1)


###############################################################################
#additional tests
#------------------------------------------------------------------------------
# 1. BAYESIAN R² 
#------------------------------------------------------------------------------

compute_bayesian_r2 <- function(fit, stan_data) {
  cat("\n=== Bayesian R² ===\n")
  
  # extract posterior samples - use rstan::extract explicitly
  posterior <- rstan::extract(fit)
  n_samples <- length(posterior$alpha)
  
  # compute fitted values for each posterior sample
  y_obs <- stan_data$mitosis_specimen
  n <- length(y_obs)
  
  # for each posterior sample, compute variance explained
  r2_samples <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    # get predicted values
    lambda_i <- exp(posterior$lambda[i, ])
    y_pred <- lambda_i
    
    # compute R²: var(fitted) / [var(fitted) + var(residuals)]
    var_fit <- var(y_pred)
    var_res <- mean((y_obs - y_pred)^2)
    r2_samples[i] <- var_fit / (var_fit + var_res)
  }
  
  r2_summary <- data.frame(
    mean = mean(r2_samples),
    sd = sd(r2_samples),
    q2.5 = quantile(r2_samples, 0.025),
    q50 = quantile(r2_samples, 0.50),
    q97.5 = quantile(r2_samples, 0.975)
  )
  
  cat("Bayesian R² Summary:\n")
  print(r2_summary)
  
  # Plot
  p <- ggplot(data.frame(r2 = r2_samples), aes(x = r2)) +
    geom_density(fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = mean(r2_samples), 
               linetype = "dashed", color = "red") +
    labs(title = "Bayesian R² Distribution",
         subtitle = sprintf("Mean R² = %.3f", mean(r2_samples)),
         x = "R²", y = "Density") +
    theme_minimal()
  
  print(p)
  
  return(list(samples = r2_samples, summary = r2_summary, plot = p))
}

#------------------------------------------------------------------------------
# 2. SIMPLE CALIBRATION PLOT 
#------------------------------------------------------------------------------

plot_simple_calibration <- function(result_part2) {
  cat("\n=== Simple Calibration Plot ===\n")
  
  pred_df <- result_part2$oos_validation$predictions
  
  # base R plot 
  plot(pred_df$observed, pred_df$predicted,
       xlab = "Observed Mitotic Count",
       ylab = "Predicted Mitotic Count (Median)",
       main = "Calibration Plot",
       pch = 19, col = as.numeric(pred_df$location),
       cex = 1.2)
  
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  
  legend("topleft", 
         legend = levels(pred_df$location),
         col = 1:4, pch = 19,
         title = "Location")
  
  # R² text
  r2_val <- cor(pred_df$observed, pred_df$predicted)^2
  text(max(pred_df$observed) * 0.1, max(pred_df$predicted) * 0.9,
       sprintf("R² = %.3f", r2_val), pos = 4)
}

#------------------------------------------------------------------------------
# 3. RESIDUAL ANALYSIS 
#------------------------------------------------------------------------------

plot_residuals <- function(result_part2) {
  cat("\n=== Residual Analysis ===\n")
  
  pred_df <- result_part2$oos_validation$predictions
  residuals <- pred_df$observed - pred_df$predicted
  
  # residual plot (like script1)
  par(mfrow = c(2, 2))
  
  # 1. Residuals vs Index
  plot(residuals, type = "h", 
       main = "Residuals vs Index",
       ylab = "Residual", xlab = "Index",
       col = as.numeric(pred_df$location))
  abline(h = 0, col = "red", lty = 2)
  
  # 2. Residuals vs Predicted
  plot(pred_df$predicted, residuals,
       main = "Residuals vs Predicted",
       xlab = "Predicted Value", ylab = "Residual",
       pch = 19, col = as.numeric(pred_df$location))
  abline(h = 0, col = "red", lty = 2)
  
  # 3. histogram of residuals
  hist(residuals, breaks = 20, 
       main = "Distribution of Residuals",
       xlab = "Residual", col = "lightblue")
  
  # 4. Q-Q plot
  qqnorm(residuals, main = "Q-Q Plot of Residuals")
  qqline(residuals, col = "red")
  
  par(mfrow = c(1, 1))
  
  # summary statistics
  cat("\nResidual Statistics:\n")
  cat(sprintf("  Mean: %.3f\n", mean(residuals)))
  cat(sprintf("  SD: %.3f\n", sd(residuals)))
  cat(sprintf("  Min: %.3f\n", min(residuals)))
  cat(sprintf("  Max: %.3f\n", max(residuals)))
}

#------------------------------------------------------------------------------
# 4. SENSITIVITY ANALYSIS 
#------------------------------------------------------------------------------

perform_sensitivity_analysis <- function(stan_code, stan_data) {
  cat("\n=== Sensitivity Analysis ===\n")
  cat("Testing robustness to prior specifications...\n")
  
  # original priors (from the paper)
  cat("\n1. Original Priors (Paper-based):\n")
  cat("  alpha ~ normal(1.75, 0.5)\n")
  cat("  delta ~ normal(1, 0.3)\n")
  cat("  epsilon ~ normal(-0.9, 0.3)\n")
  
  # fit with wider priors
  stan_code_wide <- gsub(
    "alpha ~ normal\\(1.75, 0.5\\);",
    "alpha ~ normal(1.75, 1.0);",
    stan_code
  )
  stan_code_wide <- gsub(
    "delta ~ normal\\(1, 0.3\\);",
    "delta ~ normal(1, 0.6);",
    stan_code_wide
  )
  stan_code_wide <- gsub(
    "epsilon ~ normal\\(-0.9, 0.3\\);",
    "epsilon ~ normal(-0.9, 0.6);",
    stan_code_wide
  )
  
  cat("\n2. Wider Priors (2× SD):\n")
  cat("  alpha ~ normal(1.75, 1.0)\n")
  cat("  delta ~ normal(1, 0.6)\n")
  cat("  epsilon ~ normal(-0.9, 0.6)\n")
  
  cat("\nFitting model with wider priors...\n")
  
  model_wide <- stan_model(model_code = stan_code_wide)
  fit_wide <- sampling(
    model_wide,
    data = stan_data,
    chains = 2,
    iter = 2000,
    warmup = 1000,
    refresh = 0
  )
  
  # fit with weakly informative priors
  stan_code_weak <- gsub(
    "alpha ~ normal\\(1.75, 0.5\\);",
    "alpha ~ normal(0, 10);",
    stan_code
  )
  stan_code_weak <- gsub(
    "delta ~ normal\\(1, 0.3\\);",
    "delta ~ normal(0, 10);",
    stan_code_weak
  )
  stan_code_weak <- gsub(
    "epsilon ~ normal\\(-0.9, 0.3\\);",
    "epsilon ~ normal(0, 10);",
    stan_code_weak
  )
  
  cat("\n3. Weakly Informative Priors:\n")
  cat("  alpha ~ normal(0, 10)\n")
  cat("  delta ~ normal(0, 10)\n")
  cat("  epsilon ~ normal(0, 10)\n")
  
  cat("\nFitting model with weakly informative priors...\n")
  
  model_weak <- stan_model(model_code = stan_code_weak)
  fit_weak <- sampling(
    model_weak,
    data = stan_data,
    chains = 2,
    iter = 2000,
    warmup = 1000,
    refresh = 0
  )
  
  return(list(
    fit_wide = fit_wide,
    fit_weak = fit_weak
  ))
}

compare_sensitivity_results <- function(fit_original, fit_wide, fit_weak) {
  cat("\n=== Comparing Prior Sensitivity ===\n")
  
  # extract key parameters
  params <- c("alpha", "delta", "epsilon")
  
  comparison_df <- data.frame(
    parameter = params,
    original_mean = summary(fit_original, pars = params)$summary[, "mean"],
    wide_mean = summary(fit_wide, pars = params)$summary[, "mean"],
    weak_mean = summary(fit_weak, pars = params)$summary[, "mean"]
  ) %>%
    mutate(
      diff_wide = abs(wide_mean - original_mean),
      diff_weak = abs(weak_mean - original_mean),
      max_diff = pmax(diff_wide, diff_weak)
    )
  
  print(comparison_df)
  
  cat("\nInterpretation:\n")
  if (max(comparison_df$max_diff) < 0.1) {
    cat("✓ Results are ROBUST to prior specification (max diff < 0.1)\n")
  } else if (max(comparison_df$max_diff) < 0.3) {
    cat("⚠ Results show MODERATE sensitivity (max diff 0.1-0.3)\n")
  } else {
    cat("✗ Results are SENSITIVE to priors (max diff > 0.3)\n")
  }
  
  return(comparison_df)
}

#------------------------------------------------------------------------------
# 5. K-FOLD CROSS-VALIDATION 
#------------------------------------------------------------------------------

perform_kfold_cv <- function(stan_code, full_data, K = 10) {
  cat("\n=== K-Fold Cross-Validation ===\n")
  cat(sprintf("Performing %d-fold cross-validation...\n", K))
  
  n <- nrow(full_data)
  fold_ids <- sample(rep(1:K, length.out = n))
  
  elpd_kfold <- numeric(K)
  
  for (k in 1:K) {
    cat(sprintf("Fold %d/%d...\n", k, K))
    
    # split data
    train_idx <- fold_ids != k
    test_idx <- fold_ids == k
    
    train_data <- full_data[train_idx, ]
    test_data <- full_data[test_idx, ]
    
    # prepare Stan data
    train_stan <- prepare_stan_data(train_data)
    test_stan <- prepare_stan_data(test_data)
    
    # fit model
    model <- stan_model(model_code = stan_code)
    fit <- sampling(
      model,
      data = train_stan,
      chains = 2,
      iter = 1000,
      warmup = 500,
      refresh = 0
    )
    
    # compute log predictive density on test fold
    posterior <- rstan::extract(fit)
    n_test <- test_stan$N
    
    log_pred_dens <- numeric(n_test)
    for (i in 1:n_test) {
      lambda_samples <- exp(posterior$lambda[, i])
      pred_probs <- dpois(test_stan$mitosis_specimen[i], lambda_samples)
      log_pred_dens[i] <- log(mean(pred_probs))
    }
    
    elpd_kfold[k] <- sum(log_pred_dens)
  }
  
  elpd_total <- sum(elpd_kfold)
  se_elpd <- sd(elpd_kfold) * sqrt(K)
  
  cat("\n=== K-Fold CV Results ===\n")
  cat(sprintf("ELPD: %.2f (SE: %.2f)\n", elpd_total, se_elpd))
  
  return(list(
    elpd = elpd_total,
    se = se_elpd,
    fold_elpds = elpd_kfold
  ))
}

# worflow altogether

hybrid_validation <- function(result_part1, result_part2) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  HYBRID VALIDATION: Comprehensive + Script1 Methods\n")
  cat("═══════════════════════════════════════════════════════════\n")
  
  # part 1: comprehensive validation 
  cat("\n✓ Part 1: Core Validation Complete (from result_part2)\n")
  cat("  - Convergence diagnostics\n")
  cat("  - Posterior predictive checks\n")
  cat("  - WAIC and LOO-CV\n")
  cat("  - Out-of-sample validation\n")
  cat("  - Parameter inference\n")
  
  # part 2: added methods
  cat("\n▸ Part 2: Additional Validations (from Script1)\n")
  
  # 1. bayesian R²
  cat("\n[1/4] Computing Bayesian R²...\n")
  r2_result <- compute_bayesian_r2(result_part2$fit, result_part1$train_stan)
  
  # 2. simple calibration plot
  cat("\n[2/4] Creating calibration plot...\n")
  plot_simple_calibration(result_part2)
  
  # 3. residual analysis
  cat("\n[3/4] Performing residual analysis...\n")
  plot_residuals(result_part2)
  
  # 4. sensitivity analysis
  cat("\n[4/4] Conducting sensitivity analysis...\n")
  sensitivity_fits <- perform_sensitivity_analysis(
    result_part1$stan_code,
    result_part1$train_stan
  )
  
  sensitivity_comparison <- compare_sensitivity_results(
    result_part2$fit,
    sensitivity_fits$fit_wide,
    sensitivity_fits$fit_weak
  )
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("  VALIDATION COMPLETE\n")
  cat("═══════════════════════════════════════════════════════════\n")
  
  # compile all results
  hybrid_results <- list(
    core_validation = result_part2,
    bayesian_r2 = r2_result,
    sensitivity = list(
      fits = sensitivity_fits,
      comparison = sensitivity_comparison
    )
  )
  
  return(hybrid_results)
}

################################################################################
# COMPREHENSIVE VALIDATION REPORT
################################################################################

generate_validation_report <- function(hybrid_results) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════\n")
  cat("           COMPREHENSIVE VALIDATION REPORT\n")
  cat("═══════════════════════════════════════════════════════════\n")
  
  # 1: model convergence
  cat("\n1. MODEL CONVERGENCE\n")
  cat("   ────────────────────────────────────────────────────────\n")
  
  diag <- hybrid_results$core_validation$diagnostics
  cat(sprintf("   Converged: %s\n", 
              ifelse(diag$converged, "✓ YES", "✗ NO")))
  cat(sprintf("   Max Rhat: %.4f %s\n", 
              max(diag$rhat, na.rm = TRUE),
              ifelse(max(diag$rhat, na.rm = TRUE) < 1.01, "(✓)", "(✗)")))
  cat(sprintf("   Min ESS: %.0f %s\n", 
              min(diag$ess_bulk, na.rm = TRUE),
              ifelse(min(diag$ess_bulk, na.rm = TRUE) > 400, "(✓)", "(✗)")))
  cat(sprintf("   Divergences: %d %s\n", 
              diag$divergences,
              ifelse(diag$divergences == 0, "(✓)", "(✗)")))
  
  # 2: model fit
  cat("\n2. MODEL FIT QUALITY\n")
  cat("   ────────────────────────────────────────────────────────\n")
  
  r2 <- hybrid_results$bayesian_r2$summary
  cat(sprintf("   Bayesian R²: %.3f [%.3f, %.3f]\n",
              r2$mean, r2$q2.5, r2$q97.5))
  
  ic <- hybrid_results$core_validation$information_criteria
  cat(sprintf("   WAIC: %.1f\n", 
              ic$waic$estimates["waic", "Estimate"]))
  cat(sprintf("   LOO-CV: %.1f\n", 
              ic$loo$estimates["looic", "Estimate"]))
  
  # 3: prediction accuracy
  cat("\n3. PREDICTION ACCURACY\n")
  cat("   ────────────────────────────────────────────────────────\n")
  
  oos <- hybrid_results$core_validation$oos_validation
  cat(sprintf("   MAE: %.3f\n", oos$mae))
  cat(sprintf("   RMSE: %.3f\n", oos$rmse))
  cat(sprintf("   90%% Coverage: %.3f %s\n", 
              oos$coverage_90,
              ifelse(oos$coverage_90 > 0.85 & oos$coverage_90 < 0.95, 
                     "(✓)", "(⚠)")))
  
  # 4: prior sensitivity
  cat("\n4. PRIOR SENSITIVITY\n")
  cat("   ────────────────────────────────────────────────────────\n")
  
  sens <- hybrid_results$sensitivity$comparison
  cat(sprintf("   Max parameter difference: %.3f\n", max(sens$max_diff)))
  if (max(sens$max_diff) < 0.1) {
    cat("   Assessment: ✓ ROBUST to prior specification\n")
  } else if (max(sens$max_diff) < 0.3) {
    cat("   Assessment: ⚠ MODERATELY sensitive to priors\n")
  } else {
    cat("   Assessment: ✗ SENSITIVE to prior specification\n")
  }
  
  # section 5: overall assessment
  cat("\n5. OVERALL ASSESSMENT\n")
  cat("   ────────────────────────────────────────────────────────\n")
  
  all_checks <- c(
    diag$converged,
    max(diag$rhat, na.rm = TRUE) < 1.01,
    diag$divergences == 0,
    oos$coverage_90 > 0.85 & oos$coverage_90 < 0.95,
    max(sens$max_diff) < 0.3
  )
  
  n_passed <- sum(all_checks)
  n_total <- length(all_checks)
  
  cat(sprintf("   Checks passed: %d/%d\n", n_passed, n_total))
  
  if (n_passed == n_total) {
    cat("   \n   ✓✓✓ MODEL IS READY FOR USE ✓✓✓\n")
  } else if (n_passed >= n_total - 1) {
    cat("   \n   ⚠ MODEL IS ACCEPTABLE (minor issues)\n")
  } else {
    cat("   \n   ✗ MODEL NEEDS IMPROVEMENT\n")
  }
  
  cat("\n═══════════════════════════════════════════════════════════\n\n")
}


#################
# run everything
hybrid_results <- hybrid_validation(result_part1, result_part2)

# see report
generate_validation_report(hybrid_results)

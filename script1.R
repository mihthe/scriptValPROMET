##############################
##possible validation options:

# Load libraries
library(rethinking)
library(loo)
library(bayesplot)
library(tidyverse)

# Assume 'm' as the fitted model and 'dat' as the dataset

# 1. K-Fold Cross-Validation for accuracy estimation and data generalisation
#(to be repeated for all folds)
kfold_result <- kfold(m, K = 10) #to be changed: model name, folds number
print(kfold_result)

# 2. Bayesian R² chosen to assess fit  
r2 <- bayes_R2(m)
print(r2)

# 3. Calibration Plot for predicted vs observed values plot, assess bias and calibration
predicted <- extract.samples(m)$lambda
observed <- dat$m_surg
plot(observed, apply(predicted, 2, median),
     xlab = "Observed", ylab = "Predicted Median",
     main = "Calibration Plot")
abline(0, 1, col = "red")

# 4. Residual Analysis to detect patterns or misfits
residuals <- observed - apply(predicted, 2, median)
plot(residuals, type = "h", main = "Residuals", ylab = "Residual")

# 5. Simulation-Based Calibration (simplified) to check correctness of algorithm
# Simulate prior predictive samples
prior_model <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a,
    a ~ normal(0, 1)
  ),
  data = dat,
  sample_prior = TRUE
)
prior_samples <- extract.samples(prior_model)$a
hist(prior_samples, main = "Prior Predictive Samples", xlab = "a")

# 6. Held-Out Prediction
set.seed(123)
train_index <- sample(1:nrow(dat), size = 0.8 * nrow(dat)) #training data amount changeable
train <- dat[train_index, ] 
test <- dat[-train_index, ]
#changeable for different split ratio, stratfied sampling, past vs future data splitting

held_model <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + b * Si,
    a ~ normal(0, 1),
    b ~ normal(0, 1)
  ),
  data = train
)
held_pred <- link(held_model, data = test)
held_median <- apply(held_pred, 2, median)
plot(test$m_surg, held_median,
     xlab = "Observed", ylab = "Predicted",
     main = "Held-Out Prediction")
abline(0, 1, col = "blue")

# 7. Sensitivity Analysis
# Fit model with different prior
sensitive_model <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + b * Si,
    a ~ normal(0, 5),
    b ~ normal(0, 5) #changeable priors
  ),
  data = dat
)
precis(sensitive_model)

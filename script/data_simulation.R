#libraries
library(rethinking)
library(MASS)
set.seed(2025) # to reproduce the exact figures
###########################
# functions

c.dashboard <- function (fit, warmup = FALSE, plot = TRUE, trank = TRUE) 
{
  if (class(fit) %in% c("map2stan", "ulam")) 
    fit <- fit@stanfit
  x <- rstan::get_sampler_params(fit)
  n_chains <- length(x)
  nwarmup <- fit@stan_args[[1]]$warmup
  niter <- fit@stan_args[[1]]$iter
  n_samples <- (niter - nwarmup) * n_chains
  if (warmup == FALSE) {
    for (i in 1:n_chains) {
      x[[i]] <- x[[i]][(nwarmup + 1):niter, 1:6]
    }
  }
  y <- summary(fit)
  if (n_chains > 1) {
    x_temp <- x[[1]]
    for (i in 2:n_chains) x_temp <- rbind(x_temp, x[[i]])
    x <- x_temp
  }
  else {
    x <- x[[1]]
  }
  wstart <- floor(nwarmup * 0.5)
  wend <- niter
  plot_make <- function(main, par, neff, ...) {
    ylim <- c(min(post[wstart:wend, , par]), max(post[wstart:wend, 
                                                      , par]))
    plot(NULL, xlab = "", ylab = "", type = "l", xlim = c(wstart, 
                                                          wend), ylim = ylim, ...)
    diff <- abs(ylim[1] - ylim[2])
    ylim <- ylim + c(-diff/2, diff/2)
    polygon(nwarmup * c(-1, 1, 1, -1), ylim[c(1, 1, 2, 2)], 
            col = grau(0.15), border = NA)
    neff_use <- neff
    mtext(paste("n_eff =", round(neff_use, 0)), 3, adj = 1, 
          cex = 0.9)
    mtext(main, 3, adj = 0, cex = 1)
  }
  plot_chain <- function(x, nc, ...) {
    lines(1:niter, x, col = col.alpha(rethink_palette[nc], 
                                      1), lwd = 1)
  }
  pars <- "lp__"
  post <- rstan::extract(fit, pars = pars, permuted = FALSE, inc_warmup = TRUE)
  if (plot == TRUE) {
    set_nice_margins()
    par(cex.axis = 1, cex.lab = 1.2)
    par(mfrow = c(3, 1))
    Rhat_vals <- as.numeric(round(y$summary[, 10], 2))
    plot(y$summary[, 9], Rhat_vals, xlab = "number of effective samples", 
         ylab = "Rhat", ylim = c(0.995, max(1.1, Rhat_vals, 
                                            na.rm = TRUE)))
    abline(v = 0.1 * n_samples, lty = 1, col = "red")
    abline(v = n_samples, lty = 1, col = grau())
    abline(h = 1, lty = 2)
    dens(x[, 6], adj = 0.1, xlab = "HMC energy")
    mu <- mean(x[, 6])
    sig <- sd(x[, 6])
    curve(dnorm(x, mu, sig), add = TRUE, col = rangi2, lwd = 1.5)
    if (trank == TRUE) {
      trankplot(fit, pars = "lp__", lp = TRUE, add = TRUE)
    }
    else {
      plot_make("log-probability", "lp__", y$summary["lp__", 
                                                     "n_eff"])
      for (nc in 1:n_chains) plot_chain(post[, nc, "lp__"], 
                                        nc)
    }
  }
  invisible(x)
}

##########################

# simulating N cases
N <- 1e4

# simulating U (noise - unknown causes)
U <- rnorm( n = N, mean = 0, sd = 1) 

# simulating M surg
a <-  - 3  # see Prior predictive simulation
b <-  1  # see Prior predictive simulation
log_surf <- log( 23.5 ) # surface counted in surg i.e. 5 mm^2 or 23.5 HPFs on the microscope used
mu_sur <-  log_surf + a + b * U
lambda_sur <- exp( mu_sur )
m_sur <- rpois( n = N , lambda = lambda_sur )

# simulating location
H <- 5L
L_type <- 1L:H # 1: esophageal 2: gastric 3:duodenal 4: small_int 5: large_bowel/rectum
L_prob <- c( .05 , .6 , .1 , .15 , .1 ) # probability of each site
L <- sample( x = L_type, prob = L_prob, size = N , replace = TRUE )

# simulating tumor size 
# we will use a variance covariance matrix
# we first need to simulate the elements common to the population
a <- 0 # average intercept of size standardized
b <- 0 # average slope for U of size standardized
sigma_a <- 1 # std dev intercept
sigma_b <- 1 # std slope
rho <- (0.7) # correlation between intercept and slope
# building the variance covariance matrix:
Mu <- c( a , b )
sigmas <- c( sigma_a , sigma_b ) # standard deviations
Rho <- matrix( c( 1 , rho , rho , 1 ) , nrow = 2 ) # correlation matrix
#  matrix multiply to get covariance matrix
Sigma <- diag( sigmas ) %*% Rho %*% diag( sigmas )
# and multivariate normal simulation of coefficients
corr_a_b <- mvrnorm( H , Mu , Sigma )
alpha_si <- corr_a_b[ , 1 ]
beta_si <- corr_a_b[ , 2 ]
# lastly we simulate size
mu_si <- alpha_si[L] + beta_si[L] * U
Si <- rnorm(n = N, mean = mu_si)

# simulating the counted surface of bio
# the surface of counting on the biopsy ranges form 0 to 23.5 HPFs (0 - 5 mm^2)
# it might also be a function of the location and size, being bigger location and certain size
# easier to be biopsied

# we will use a variance covariance matrix
# we first need to simulate the elements common to the population
a <- 0 # average intercept of inv_log surface standardized
b <- 0 # average slope for size of inv_log surface standardized
sigma_a <- 1 # std dev intercept
sigma_b <- 1 # std slope
rho <- (0.7) # correlation between intercept and slope
# building the variance covariance matrix:
Mu <- c( a , b )
sigmas <- c( sigma_a , sigma_b ) # standard deviations
Rho <- matrix( c( 1 , rho , rho , 1 ) , nrow = 2 ) # correlation matrix
#  matrix multiply to get covariance matrix
Sigma <- diag( sigmas ) %*% Rho %*% diag( sigmas )
# and multivariate normal simulation of coefficients
corr_a_b <- mvrnorm( H , Mu , Sigma )
alpha_surf <- corr_a_b[ , 1 ]
beta_surf <- corr_a_b[ , 2 ]
# lastly we simulate size
mu_surf <- 2 + alpha_surf[L] + beta_surf[L] * Si
logit_Su <- rnorm(n = N, mean = mu_surf, sd = 2)
Su <- inv_logit(logit_Su)
Su <- Su * 23.5 

# now we can simulate M bio
a <- -3 # lambd_sur and log sur add + 4   
b <- 1
log_surf <- log( Su )
mu_bio <-  log_surf + a + b * mu_sur
lambda_bio <- exp( mu_bio )
m_bio <- rpois( n = N , lambda = lambda_bio )

# lastly let simulate therapy and response to it
# therapy in GIST is done in high risk patients
# risk is based on size site and mitosis
# we will create cutoffs
# for size in each location
H <- length(unique(L))
si_L_cutoff <- vector(length = H)
for(i in 1:H) si_L_cutoff[i] <- quantile(Si[L==i], 0.5) 
# and for mitotic count on biopsy
m_bio_cutoff <- qpois( 0.5, mean(m_bio)) 
# next we calculate the risk class and assume that all the
# high risk (R = 1) will undergo neoadjuvant therapy
R <- m_bio > m_bio_cutoff & Si > si_L_cutoff[L]

# in our simulation response to therapy will be inversely
# associated with biological aggressiveness
# more aggressive cases will respond less
p_resp_Rx <- 1 - inv_logit( U - 1.5 )
resp_Rx <- vector(mode = 'integer',  length = N)
for(i in which(R==TRUE)) resp_Rx[i] <- rbinom(n = 1, size = 1, p_resp_Rx[i])
for(i in which(resp_Rx==1)) m_sur[i] <- rpois(n = 1, lambda = 0.1)

# build the db with N cases
sim <- data.frame( m_bio , m_sur, Su , Si , L , U, resp_Rx)

# creating a slim db for the analysis
slim <- sample(x = 1:nrow(sim) , size = 200) #index to sample
train <- slim[1:round(length(slim)/2)]
validation <- slim[round(length(slim)/2):length(slim)]
slim_sim <- sim[train, ]
pairs(slim_sim)

##### Data analysis 

dat <- list(
  #observed data
  Si = slim_sim$Si,             # size
  Su = standardize(slim_sim$Su),      # Biopsy surface
  m_bio = standardize(slim_sim$m_bio),         #mitosis on biopsy
  m_surg = slim_sim$m_sur,       #mitosis on surgery
  L = slim_sim$L,        # possible location of GIST "gastric","duodenum","ileum","colon","rectum", "esophagus"
  n = as.integer(1 - slim_sim$resp_Rx),
  y = slim_sim$resp_Rx
)

# fitting the model
# m <- ulam(
#   alist(
#     m_surg ~ dpois(lambda),
#     log(lambda) <- a + b[L] * Si + g * m_bio * n  + d * m_bio * y + e[L] * Su,
#     b[L] ~ normal( b_bar , sigmab),
#     e[L] ~ normal( e_bar , sigmae),
#     a ~ normal( -1 , 0.2),
#     c(b_bar, g , d , e_bar ) ~ normal( -1 , 0.2 ),
#     c(sigmab, sigmae) ~ dexp( 2 )
#   ), data=dat , chains=4 , cores=4 , iter=4000 , log_lik = FALSE)
# 2% of divergent transitions therefore we moved to NC
# fitting the model nc
m <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + sigmab * z_b[L] * Si + g * m_bio * n  + d * m_bio * y + sigmae * z_e[L] * Su,
    z_b[L] ~ normal( -1 , 0.5),
    z_e[L] ~ normal( -1 , 0.5),
    a ~ normal( -1 , 0.2),
    c( g , d  ) ~ normal( -1 , 0.2 ),
    c(sigmab, sigmae) ~ dexp( 2 )
  ), data=dat , chains=4 , cores=4 , iter=8000 , log_lik = FALSE)

# model diagnostics
jpeg(paste0("output/figures/","dashboardSIM",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
c.dashboard(m)
dev.off()

jpeg(paste0("output/figures/","trankplotSIM",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
trankplot(m, n_cols = 4)
dev.off()

# Posterior probability
post <- extract.samples(m) # extract the samples from the model's fit
precis(m,2) # makes a table with model coefficients

# plot the model coefficients
jpeg(paste0("output/figures/","Model coefficients sim",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
plot(NULL, xlim = c(- 1.3, 1.1), ylim = c(0,11), xlab = '', 
     ylab= 'Density', main = 'Model coefficients')
abline(v = 0, col = scales::alpha(1,0.4))
for(i in 2:5) dens(post[[i]], add = TRUE, lwd= 3, col = i-1)
for(i in 1:5) dens(post$z_b[,i], add = TRUE, lwd= 3, col = 4+i, lty = 2)
legend('topleft',lwd = 3, col = 1:8, lty = c(1,1,1,1,2,2,2,2,2),
       legend = c('bio_surface slope', 'intercept', 'm_bio_NAC_resp slope', 'm_bio slope', 
                  paste0('size slope [',
                        c('esophagus', 'stomach', 'duodenum','small bowel','colonrectum'),']')))
dev.off()

# Posterior predictive check
# the variables were standardized* therefore to make a retrodiction
# we need to transform them in the standardized 
# *: divided by the standard deviation and subtracted by the mean
mubio <- attr(dat$m_bio, 'scaled:center')
sdbio <- attr(dat$m_bio, 'scaled:scale')
muSu <- attr(dat$Su, 'scaled:center')
sdSu <- attr(dat$Su, 'scaled:scale')

# the 'link' function will return a lambda distribution given the size, the site, the surface, 
# and the mitotic count on biopsy
# m_link <- function( Si , m_bio, Su , L ) {
#   m_bio <- (m_bio - mubio) / sdbio
#   Su <- (Su - muSu) / sdSu
#   mu <- with( post ,{
#     a + b[,L] * Si + g * m_bio + e[,L] * Su })
#   lambda <- exp( mu )
#   lambda
# }


m_link <- function( Si , m_bio, Su , L ) {
  m_bio <- (m_bio - mubio) / sdbio
  Su <- (Su - muSu) / sdSu
  mu <- with( post ,{
    a + sigmab * z_b[ ,L] * Si + g * m_bio + sigmae * z_e[ ,L] * Su })
  lambda <- exp( mu )
  lambda
}

#as an example
lambda <- m_link( Si =  5 , m_bio = 7, Su = 23.5, L = 2)
median(lambda)
HPDI(lambda)

# therefore to make the posterior predictive check we used
# our synthetic dataset 
lambda <- mapply(m_link , Si = slim_sim$Si, m_bio = slim_sim$m_bio , 
                 Su = slim_sim$Su , L = slim_sim$L )
lmed <- apply( lambda , 2 , median ) #
lci <- apply( lambda , 2 , HPDI )
j <- order(lmed)
# 
# jpeg(paste0("output/figures/","PosPredCheck_sim_only_lambda",".jpg"), 
#      units = "in", 
#      width = 15, height = 7, res = 300)
# par(mfrow=c(1,1))
# plot(NULL, xlim = c(1,length(slim)), ylim = c(log(0.1),log(70)), bty = 'n', xaxt = 'n', 
#      yaxt = "n", ylab = 'lambda: expected value of mitotic count', 
#      xlab = 'Cases', main = 'Posterior predictive check')
# abline(h = log(c(0.1,1,2,3,5,10,20, 50, 100)), col = scales::alpha(1,0.3))
# points(1:length(slim), log(lambda_sur[slim][j]), col = 6, pch = 4, lwd= 2)
# segments(x0 = 1:length(slim), y0 = log(lci[1,j]+ 0.1), y1 = log(lci[2,j]+ 0.1), lwd = 2)
# points(1:length(slim), log(lmed[j] + 0.1), pch = 16)
# axis(2, at = log(c(0.1,1,2,3,5,10,20, 50, 100)), labels = c(0.1,1,2,3,5,10,20,50, 100), las = 2)
# axis(1, at = 1:length(slim), labels = 1:length(slim))
# legend('topleft', legend = c('True lambda surgery','Predicted lambda surgery'), 
#        pch = c(4, 16), lwd = 2, lty = c(0 ,1), col = c(6, 1))
# dev.off()
# 
# jpeg(paste0("output/figures/","PosPredCheck_sim",".jpg"), 
#      units = "in", 
#      width = 15, height = 7, res = 300)
# par(mfrow = c(1,1))
# plot(NULL, xlim = c(1,length(slim)), ylim = c(log(0.1),log(70)), bty = 'n', xaxt = 'n', 
#      yaxt = "n", ylab = 'Mitotic count', xlab = 'Cases', main = 'Posterior predictive check')
# abline(v = 1:length(slim), lty = 1, col = scales::alpha(dat$L[j], 0.2), lwd = 10)
# abline(h = log(5+ 0.1))
# points(1:length(slim), log(slim_sim$m_bio[j] + 0.1), col = 4, pch = 4, lwd= 2, cex = 0.3 + slim_sim$Su/23.5)
# points(1:length(slim), log(dat$m_surg[j] + 0.1), col = 2 + dat$n, lwd= 3, cex = 0.3 + inv_logit(dat$Si) )
# segments(x0 = 1:length(slim), y0 = log(lci[1,j]+ 0.1), y1 = log(lci[2,j]+ 0.1), lwd = 2)
# points(1:length(slim), log(lmed[j] + 0.1), pch = 16)
# axis(2, at = log(c(0.1,1,2,3,5,10,20, 50, 100)), labels = c(0,1,2,3,5,10,20,50, 100), las = 2)
# axis(1, at = 1:length(slim), labels = 1:length(slim))
# legend('topleft', legend = c('Biopsy','Surgery','Surgery NAC response','Prediction'), 
#        pch = c(4, 1, 1, 16), lwd = 2, lty = c(0, 0 , 0 ,1), col = c(4, 3, 2, 1))
# dev.off()


jpeg(paste0("output/figures/","PosPredCheck_sim_full",".jpg"), 
     units = "in", 
     width = 15, height = 15, res = 300)
par(mfrow=c(2,1))
plot(NULL, xlim = c(1,length(train)), ylim = c(log(0.1),log(70)), bty = 'n', xaxt = 'n', 
     yaxt = "n", ylab = 'lambda: expected value of mitotic count', 
     xlab = 'Cases', main = 'Posterior predictive check')
abline(h = log(c(0.1,1,2,3,5,10,20, 50, 100)), col = scales::alpha(1,0.3))
points(1:length(train), log(lambda_sur[slim][j]), col = 6, pch = 5, lwd= 2)
segments(x0 = 1:length(train), y0 = log(lci[1,j]+ 0.1), y1 = log(lci[2,j]+ 0.1), lwd = 2)
points(1:length(train), log(lmed[j] + 0.1), pch = 16)
axis(2, at = log(c(0.1,1,2,3,5,10,20, 50, 100)), labels = c(0.1,1,2,3,5,10,20,50, 100), las = 2)
axis(1, at = 1:length(train), labels = 1:length(train))
legend('topleft', legend = c('True lambda surgery','Predicted lambda surgery'), 
       pch = c(5, 16), lwd = 2, lty = c(0 ,1), col = c(6, 1))

plot(NULL, xlim = c(1,length(train)), ylim = c(log(0.1),log(70)), bty = 'n', xaxt = 'n', 
     yaxt = "n", ylab = 'Mitotic count', xlab = 'Cases', main = 'Posterior predictive check')
abline(v = 1:length(train), lty = 1, col = scales::alpha(dat$L[j], 0.2), lwd = 10)
abline(h = log(5+ 0.1))
points(1:length(train), log(slim_sim$m_bio[j] + 0.1), col = 4, pch = 4, lwd= 2, cex = 0.3 + slim_sim$Su/23.5)
points(1:length(train), log(dat$m_surg[j] + 0.1), col = 2 + dat$n, lwd= 3, cex = 0.3 + inv_logit(dat$Si) )
segments(x0 = 1:length(train), y0 = log(lci[1,j]+ 0.1), y1 = log(lci[2,j]+ 0.1), lwd = 2)
points(1:length(train), log(lmed[j] + 0.1), pch = 16)
axis(2, at = log(c(0.1,1,2,3,5,10,20, 50, 100)), labels = c(0,1,2,3,5,10,20,50, 100), las = 2)
axis(1, at = 1:length(train), labels = 1:length(train))
legend('topleft', legend = c('Biopsy','Surgery','Surgery NAC response','Prediction'), 
       pch = c(4, 1, 1, 16), lwd = 2, lty = c(0, 0 , 0 ,1), col = c(4, 3, 2, 1))
dev.off()


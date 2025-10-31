# libraries
library(tidyverse)
library(rethinking)
set.seed(2022)
###########################
# functions
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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

###########################
# importing datasets
db <- read.csv("db_upload.csv")
# 
# where Si = Size, Su = Surface, m_bio = mitotic count on biopsy, 
# m_sur = mitotic count on surgery. Si, Su, m_bio were standardize()

###########################

dat <- as.list(db)


### the model fit
# to chose among the best model we compared several alternative models using
# WAIC and PSIS-LOO

h_size_surf <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + b[L] * Si + g * m_bio * n  + d * m_bio * y + e[L] * Su,
    b[L] ~ normal( b_bar , sigmab),
    e[L] ~ normal( e_bar , sigmae),
    a ~ normal( -1 , 0.2),
    c(b_bar, g , d , e_bar ) ~ normal( -1 , 0.2 ),
    c(sigmab, sigmae) ~ dexp( 2 )
  ), data=dat , chains=4 , cores=4 , iter=4000 , log_lik = TRUE)

h_size <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + b[L] * Si + g * m_bio * n  + d * m_bio * y + e * Su,
    b[L] ~ normal( b_bar , sigma),
    a ~ normal( -1 , 0.2),
    c( b_bar , g , d , e ) ~ normal( -1 , 0.2 ),
    sigma ~ dexp( 2 )
  ), data=dat , chains=4 , cores=4 , iter=4000 , log_lik = TRUE)

h_size_surf_m_bio <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + b[L] * Si + g[L] * m_bio * n  + d * m_bio * y + e[L] * Su,
    b[L] ~ normal( b_bar , sigmab),
    g[L] ~ normal( g_bar , sigmag),
    e[L] ~ normal( e_bar , sigmae),
    a ~ normal( -1 , 0.2),
    c(b_bar, g , d , e_bar, g_bar ) ~ normal( -1 , 0.2 ),
    c(sigmab, sigmae, sigmag) ~ dexp( 2 )
  ), data=dat , chains=4 , cores=4 , iter=4000 , log_lik = TRUE)

no_size  <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a +  g * m_bio * n  + d * m_bio * y + e * Su,
    a ~ normal( -1 , 0.2),
    c(  g , d , e ) ~ normal( -1 , 0.2 )
  ), data=dat , chains=4 , cores=4 , iter=4000 , log_lik = TRUE)

h_size_no_surf <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + b[L] * Si + g * m_bio * n  + d * m_bio * y,
    b[L] ~ normal( b_bar , sigma),
    a ~ normal( -1 , 0.2),
    c( b_bar , g , d ) ~ normal( -1 , 0.2 ),
    sigma ~ dexp( 2 )
  ), data=dat , chains=4 , cores=4 , iter=4000 , log_lik = TRUE)

##### model comparison

jpeg(paste0("output/figures/","WAIC comparison",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
plot(compare(h_size_surf, h_size, h_size_surf_m_bio, 
             no_size, h_size_no_surf, func = "WAIC"))
dev.off()

jpeg(paste0("output/figures/","LOO comparison",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
plot(compare(h_size_surf, h_size, h_size_surf_m_bio, 
             no_size, h_size_no_surf, func = "LOO"))
dev.off()
########## Prior predictive simulation

m <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + b[L] * Si + g * m_bio * n  + d * m_bio * y + e[L] * Su,
    b[L] ~ normal( b_bar , sigmab),
    e[L] ~ normal( e_bar , sigmae),
    a ~ normal( -1 , 0.2),
    c(b_bar, g , d , e_bar ) ~ normal( -1 , 0.2 ),
    c(sigmab, sigmae) ~ dexp( 2 )
  ), data=dat , chains=4 , cores=4 , iter=4000 , sample_prior = TRUE)


post <- extract.samples(m)

jpeg(paste0("output/figures/","Model coefficients gastric PPS",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
plot(NULL, xlim = c(-2, 1.5), ylim = c(0,2.3), xlab = '', 
     ylab= 'Density', main = 'Model coefficients', sub= '*: displayed only for Gastric GIST')
abline(v = 0, col = scales::alpha(1,0.4))
for(i in 1:2) dens(post[[i]][,4], add = TRUE, lwd= 3, col = i)
for(i in 5:6) dens(post[[i]], add = TRUE, lwd= 3, col = i-2)
legend('topright',lwd = 3, col = 1:4, legend = c('tumor size *', 'biopsy surface*', 
                                                 'biosy count (if NAC_resp)', 'biopsy count'))
dev.off()

m_link <- function( Si , m_bio, Su , L ) {
  Si <- (Si - 62.175) / 48.07636
  m_bio <- (m_bio - 3.6625) / 8.066094
  Su <- (Su -14.51) / 9.335991
  mu <- with( post ,{
    a + b[,L] * Si + g * m_bio + e[,L] * Su })
  lambda <- exp( mu )
  lambda
}
lambda <- m_link( Si =  100 , m_bio = 7, Su = 23.5, L = 2)
median(lambda)
HPDI(lambda)

lambda <- mapply(m_link , Si = db$size, m_bio = db$m_bio , Su = db$Su , L = dat$L )
lmu <- apply( lambda , 2 , median )
lci <- apply( lambda , 2 , HPDI )
j <- order(lmu)

jpeg(paste0("output/figures/","PriorPredCheck",".jpg"), 
     units = "in", 
     width = 12, height = 7, res = 300)
par(mfrow = c(1,1))
plot(NULL, xlim = c(1,80), ylim = c(log(0.1),log(2000)), 
     bty = 'n', xaxt = 'n',  yaxt = "n", ylab = 'Mitotic count', 
     xlab = 'Cases', main = 'Prior predictive simulation')
abline(h = log(5+ 0.1))
segments(x0 = 1:80, y0 = log(lci[1,j]+ 0.1), 
         y1 = log(lci[2,j]+ 0.1), lwd = 2)
points(1:80, log(lmu[j] + 0.1), pch = 16)
axis(2, at = log(c(0.1,1,2,3,5,10,20, 50, 100)), 
     labels = c(0,1,2,3,5,10,20,50, 100), las = 2)
axis(1, at = 1:80, labels = 1:80)
dev.off()


########## re-computing w/ out log_lik for trankplot saving

m <- ulam(
  alist(
    m_surg ~ dpois(lambda),
    log(lambda) <- a + b[L] * Si + g * m_bio * n  + d * m_bio * y + e[L] * Su,
    b[L] ~ normal( b_bar , sigmab),
    e[L] ~ normal( e_bar , sigmae),
    a ~ normal( -1 , 0.2),
    c(b_bar, g , d , e_bar ) ~ normal( -1 , 0.2 ),
    c(sigmab, sigmae) ~ dexp( 2 )
  ), data=dat , chains=4 , cores=4 , iter=4000 , log_lik = FALSE)

jpeg(paste0("output/figures/","dashboard",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
c.dashboard(m)
dev.off()

jpeg(paste0("output/figures/","trankplot",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
trankplot(m, n_cols = 4)
dev.off()

#the following code export it as a rmarkdown table
coeff_table <- precis(m,2)
table <- matrix(unlist(coeff_table@.Data), nrow = 15) |> round(2)
table[,5] <- table[,5] |> round()
coeff_table@row.names
# in the paper we ordered the labels so we need to change epsilon became gamma, gamma became delta, and delta became epsilon
newRowName <-  c("$\\beta_{L[1]}$","$\\beta_{L[2]}$","$\\beta_{L[3]}$","$\\beta_{L[4]}$",
                 "$\\gamma_{L[1]}$","$\\gamma_{L[2]}$", "$\\gamma_{L[3]}$","$\\gamma_{L[4]}$",
                 "$\\alpha$","$\\bar{\\gamma}$" , "$\\epsilon$","$\\delta$", "$\\bar{\\beta}$", 
                 "$\\sigma_{\\gamma}$","$\\sigma_{\\beta}$")
coeff_table@names
newColName <- c("mean", "SD", "5.5%", "94.5%", "N eff.",  "$\\hat{R}$")
row.names(table) <-  newRowName
colnames(table) <-  newColName
knitr::kable(table, caption = '(\\#tab:coefValue) Coefficients')


# extract the posterior probability 
post <- extract.samples(m)
# # the following command save the posterior probability density of the coefficients
# saveRDS(post, '~/Dropbox/R/R_Projects/Surprise_GIST/mit_Cal/post.rds')


jpeg(paste0("output/figures/","Model coefficients gastric",".jpg"), 
     units = "in", 
     width = 7, height = 6, res = 300)
plot(NULL, xlim = c(-0.35, 1.5), ylim = c(0,12), xlab = '', 
     ylab= 'Density', main = 'Model coefficients', sub= '*: displayed only for Gastric GIST')
abline(v = 0, col = scales::alpha(1,0.4))
for(i in 1:2) dens(post[[i]][,4], add = TRUE, lwd= 3, col = i)
for(i in 5:6) dens(post[[i]], add = TRUE, lwd= 3, col = i-2)
legend('topright',lwd = 3, col = 1:4, legend = c('tumor size *', 'biopsy surface*', 
                                                 'biosy count (if NAC_resp)', 'biopsy count'))
dev.off()

m_link <- function( Si , m_bio, Su , L ) {
  Si <- (Si - 62.175) / 48.07636
  m_bio <- (m_bio - 3.6625) / 8.066094
  Su <- (Su -14.51) / 9.335991
  mu <- with( post ,{
    a + b[,L] * Si + g * m_bio + e[,L] * Su })
  lambda <- exp( mu )
  lambda
}
lambda <- m_link( Si =  100 , m_bio = 7, Su = 23.5, L = 2)
median(lambda)
HPDI(lambda)

lambda <- mapply(m_link , Si = db$size, m_bio = db$m_bio , Su = db$Su , L = dat$L )
lmed <- apply( lambda , 2 , median )
lci <- apply( lambda , 2 , HPDI )
j <- order(lmed)

jpeg(paste0("output/figures/","PosPredCheck",".jpg"), 
     units = "in", 
     width = 12, height = 7, res = 300)
par(mfrow = c(1,1))
plot(NULL, xlim = c(1,80), ylim = c(log(0.1),log(2000)), 
     bty = 'n', xaxt = 'n',  yaxt = "n", ylab = 'Mitotic count', 
     xlab = 'Cases', main = 'Posterior predictive check')
abline(v = 1:80, lty = 1, col = scales::alpha(dat$L[j], 0.2), lwd = 10)
abline(h = log(5+ 0.1))
points(1:80, log(db$m_bio[j] + 0.1), col = 4, pch = 4, 
       lwd= 2, cex = 0.3 + db$Su/23.5)
points(1:80, log(dat$m_surg[j] + 0.1), col = 2 + dat$n[j], 
       lwd= 3, cex = 0.3 + inv_logit(dat$Si) )
segments(x0 = 1:80, y0 = log(lci[1,j]+ 0.1), 
         y1 = log(lci[2,j]+ 0.1), lwd = 2)
points(1:80, log(lmed[j] + 0.1), pch = 16)
axis(2, at = log(c(0.1,1,2,3,5,10,20, 50, 100)), 
     labels = c(0,1,2,3,5,10,20,50, 100), las = 2)
axis(1, at = 1:80, labels = 1:80)
legend('topleft', legend = c('Biopsy','Surgery','Surgery NAC response',
                             'Prediction', "Stomach", "Duodenum", 
                             'Small-bowel', 'Colon-rectum'), 
       pch = c(4, 1, 1, 16, 22, 22, 22, 22), 
       lwd = c(2, 2, 2, 2, 0, 0, 0, 0), 
       lty = c(0, 0 , 0 ,1, 0, 0, 0, 0), 
       col = c(4, 3, 2, 1, 1, 1, 1, 1),
       pt.cex = c(1, 1, 1, 1, 2.5, 2.5, 2.5, 2.5),
       pt.bg = c( 0, 0, 0, 0, 
                  scales::alpha(c(4, 3, 2, 1), 0.3)), ncol = 2)
dev.off()
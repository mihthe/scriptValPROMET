set.seed(260109)

N <- 100 # number of pts
l_s <- 1e4 # length of simulation


Y_t <- rpois(n = N, lambda = 2) # simulate the count on the surgical specimen
p_e <- rexp(n = N, rate = 5) # simulate the error of the biopsy count
Y_p <- matrix(nrow = N, ncol = l_s) # initialize the matrix for the simulated posterior

#populate the matrix of the simulated posterior
for(i in 1:N) {
  Y_p[i,] <- rpois(n = l_s, lambda = Y_t[i] + p_e)
}

resid <- matrix(nrow = N, ncol = l_s) # initialize the matrix for the residual 

#populate the matrix
for(i in 1:N) {
  resid[i,] <- Y_t[i] - Y_p[i,]
}

hist(resid, freq = FALSE, 
     main = "Residuals", 
     xlab = "Error")
####################
### Residuals CI ###
####################
rethinking::HPDI(c(resid))


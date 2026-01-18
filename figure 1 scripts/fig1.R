### for fig1. 
# this, step 4

# library
library(scales)

# lambda from step 3
# slim_sim (subset of simulated dataset) from step 1
# m_sur from step 1

N <- nrow(slim_sim)  # total number of patients (rows of your dataset)

# lambda for slim_sim dataset of 100 patients
lambda <- mapply(m_link, Si = slim_sim$Si, m_bio = slim_sim$m_bio, 
                 Su = slim_sim$Su, L = slim_sim$L)
# patients by median predicted lambda 
lambda_median <- apply(lambda, 2, median)  
idx <- order(lambda_median)  

# y-axis limits based on true mitotic counts (m_sur) and predicted lambda
y_min <- min(c(slim_sim$m_sur, 0))
y_max <- max(c(slim_sim$m_sur, apply(lambda, 2, quantile, probs = 0.95)))

# plot as JPEG image
jpeg(paste0("output/figures/", "fig1.jpg"), 
     units = "in", 
     width = 14, height = 8, res = 300)

par(mai = c(0.9, 0.9, 0.9, 0.4))

plot(NULL, 
     xlim = c(0.5, N + 0.5),  
     ylim = c(y_min, y_max),  
     xlab = 'Patient from simulation subset', 
     ylab = 'Mitotic count',
     main = 'Figure 1: Predictive Check',
     bty = 'l')

# reference line at mitotic count = 5
abline(h = 5, lty = 2, col = "gray40", lwd = 2)

# scaling factor for violin width ì
s_factor <- 0.4  

# through each patient (in ordered sequence)
for(i in 1:N) {
  patient_idx <- idx[i]  # the actual patient index
  
  # posterior distribution for this patient (lambda)
  patient_lambda <- lambda[, patient_idx]
  
  # density of the posterior distribution
  dens_obj <- density(patient_lambda, adjust = 1.2)
  x_dens <- dens_obj$y  # density values (width of violin)
  y_dens <- dens_obj$x  # lambda values (height of violin)
  
  # normalize density for plotting (avoiding overlap)
  x_dens_norm <- x_dens / max(x_dens) * s_factor
  
  # right half of violin
  polygon(i + x_dens_norm, y_dens, 
          col = scales::alpha("skyblue", 0.5), 
          border = scales::alpha("steelblue", 0.7), 
          lwd = 0.3)
  
  # left half of violin (mirror)
  polygon(i - x_dens_norm, y_dens, 
          col = scales::alpha("skyblue", 0.5), 
          border = scales::alpha("steelblue", 0.7), 
          lwd = 0.3)
  
  # mitotic count on surgery specimen (m_sur)
  points(i, slim_sim$m_sur[patient_idx], 
         pch = 4, col = "blue", lwd = 2, cex = 1.2)  # blue crosses for m_sur
  
  # median prediction as a horizontal line
  segments(i - 0.35, lambda_median[patient_idx], 
           i + 0.35, lambda_median[patient_idx], 
           lwd = 2, col = "black")
}

# legend
legend('topleft', 
       pch = c(4, 15),  # Blue crosses for m_sur
       col = c("blue", "lightblue"),
       pt.cex = c(1.2, 2),
       legend = c("Surgical specimen (m_sur)", 
                  "Posterior Predictive Distribution (lambda)"),
       bty = "n")

# text
mtext("Violin shapes show posterior predicted distribution", 
      side = 3, line = 0.5, cex = 0.8, adj = 0)

dev.off()
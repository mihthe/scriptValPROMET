sim_sur_l <- rpois(n = length(lambda), lambda = lambda)
sim_sur <- matrix(sim_sur_l, nrow = nrow(lambda), ncol = ncol(lambda))

# raw and absolute differences
error <- sim_sur_l - db$m_sur
error <- matrix(error, nrow = nrow(lambda))
db$error <- colMeans(error)

jpeg(paste0("output/figures/", "Posterior_Difference_validation.jpg"), 
     units = "in", 
     width = 5, height = 5, res = 400)
plot(density(error), xlim = c(-30,30), lwd = 3, 
main = 'Signed Prediction Error on Validation set', 
xlab = "Posterior Probability of Signed Difference (sim_surg - m_surg)", adj=0.5)
abline(v=0,lty=2)
legend("topleft", legend = "model underestimate", bty="n")
legend("topright", legend = "model overestimate", bty="n")
dev.off()
# post difference 
idx <- order(db$error)

jpeg(paste0("output/figures/", "Posterior_Difference_casewise_validation.jpg"), 
     units = "in", 
     width = 10, height = 5, res = 400)
plot(NULL, ylim=c(-30,30), xlim = c(0.5,34.5),
     xlab = 'Cases', ylab = 'Density', xaxt = 'null', 
     main = 'Posterior Probability of Signed Difference (sim_surg - m_surg)')
abline(h = 0)
s_factor <- 0.3 # scaling factor for graphics
for(i in 1:34) {
  y <- density(error[,idx[i]])$x
  x <- density(error[,idx[i]])$y
  polygon(i + x/s_factor, y, col = scales::alpha(db$Site[idx[i]] ,0.6), border = FALSE)
  lines(i + x/s_factor, y, lwd = 1)
  polygon(i - x/s_factor, y, col = scales::alpha(db$Site[idx[i]],0.6), lwd = 2, border = FALSE)
  lines(i - x/s_factor, y, lwd = 1)
}
dev.off()

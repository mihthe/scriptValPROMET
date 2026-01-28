post <- readRDS("PATH/post.rds")

str(post)

link <- function( Si , m_bio , Su , L ) {
 
  # transform the measure to the 'standardized' version (mean = 0 and sd = 1)
  # i.e. subtracted by the population mean and divided by the standard deviation
  
  Si <- (Si - 62.175) / 48.07636          
  m_bio <- (m_bio - 3.6625) / 8.066094
  Su <- (Su -14.51) / 9.335991
  
  # construct a vector of the log(lambda)
  mu <- with( post ,{
    a + b[,L] * Si + g * m_bio + e[,L] * Su })
  
  # convert it in to positve value 
  lambda <- exp( mu )
  lambda
} 


### validating on the remaining
# our synthetic dataset 
lambda <- mapply(link , Si = db$Size , m_bio = db$m_bio , 
                 Su = db$Su , L = db$Site)
lmed <- apply( lambda , 2 , median ) #
lci <- apply( lambda , 2 , rethinking::HPDI )
j <- order(lmed)
# 

jpeg(paste0("output/figures/","Figure1_validation",".jpg"), 
     units = "in", 
     width = 7, height = 5, res = 300)
plot(NULL, xlim = c(1,length(db[,1])), ylim = c(log(0.1),log(150)), bty = 'n', xaxt = 'n', 
     yaxt = "n", ylab = 'Mitotic count', xlab = 'Cases', main = 'Predictive check')
abline(v = 1:length(db[,1]), lty = 1, col = scales::alpha(db$Site[j], 0.2), lwd = 10)
abline(h = log(5+ 0.1))
points(1:length(db[,1]), log(db$m_bio[j] + 0.1), col = 4, pch = 4, lwd= 2, cex = 0.3 + db$Su/23.5)
points(1:length(db[,1]), log(db$m_sur[j] + 0.1), col = 2 , lwd= 3, cex = 0.3 + db$Su/23.5 )
segments(x0 = 1:length(db[,1]), y0 = log(lci[1,j]+ 0.1), y1 = log(lci[2,j]+ 0.1), lwd = 2)
points(1:length(db[,1]), log(lmed[j] + 0.1), pch = 16)
axis(2, at = log(c(0.1,1,2,3,5,10,20, 50, 100)), labels = c(0,1,2,3,5,10,20,50, 100), las = 2)
axis(1, at = 1:length(db[,1]), labels = 1:length(db[,1]))
legend('topleft', legend = c('Biopsy','Surgery','Prediction'), 
       pch = c(4, 1, 16), lwd = 2, lty = c(0, 0 ,1), col = c(4, 2, 1))
legend('bottomright', legend = c("Stomach", "Duodenum", 'Small-bowel', 'Colon-rectum'),
pch = 22, pt.bg = scales::alpha(c(4,2,3,1), 0.5))       
dev.off()

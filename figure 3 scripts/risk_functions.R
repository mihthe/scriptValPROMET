# run data simulation file

# save datasets from simulation file
write.csv(sim, "sim.csv", row.names = FALSE)
write.csv(slim_sim, "slim_sim.csv", row.names = FALSE)

# run data analysis file 

# risk function from app file 
risk_str <- function(size, mic, site) {
  
  size <- size / 10
  
  if (site == 4L) {
    
    if (mic <= 5) {
      if (size <= 2) {
        return("none")
      } else if (size <= 5) {
        return("very low")
      } else if (size <= 10) {
        return("low") 
      } else {
        return("moderate")
      }
      
    } else { # mic > 5
      if (size <= 2) {
        return("none")
      } else if (size <= 5) {
        return("moderate")
      } else if (size <= 10) {
        return("high") 
      } else {
        return("high")
      }
    }
    
  } else { 
    
    if (mic <= 5) {
      if (size <= 2) {
        return("none")
      } else if (size <= 5) {
        return("low")
      } else if (size <= 10) {
        return("moderate") 
      } else {
        return("high")
      }
      
    } else { # mic > 5
      if (size <= 2) {
        return("high")
      } else if (size <= 5) {
        return("high")
      } else if (size <= 10) {
        return("high") 
      } else {
        return("high")
      }
    }
  }
}

# Vectorized version for applying to entire dataframe or vector
risk_str_vectorized <- function(size, mic, site) {
  
  # handle single values and vectors
  n <- max(length(size), length(mic), length(site))
  
  # initialize output vector
  risk <- character(n)
  
  # risk function to each case
  for (i in 1:n) {
    risk[i] <- risk_str(
      size = size[i], 
      mic = mic[i], 
      site = site[i]
    )
  }
  
  return(risk)
}

##################################
# to use risk functions elsewhere
##################################

# at the top of R script:
source("risk_functions.R")

# if a single patient:
risk_single <- risk_str(size = 45, mic = 7, site = 4L) # arbitrary values

# if multiple patients (where db is your dataset):
db$risk_biopsy <- risk_str_vectorized(
  size = db$Si, 
  mic = db$m_bio, 
  site = db$L
)
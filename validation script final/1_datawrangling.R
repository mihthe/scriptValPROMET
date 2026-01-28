# this script is for data wrangling
# importing the data
db <- read.csv("~/Downloads/db.csv")
str(db)
na_idx <- db[,] == ""
db[na_idx] <- NA
unique(db$Insert_Patient_ID)
mean((db$Is_it_a_biopsy. == "yes") == (!is.na(db$biopsy_N)))
mean((db$Is_it_a_biopsy. == "no") == (!is.na(db$surgical_N)))
naj_na <- is.na(db$neoadjuvant_therapy.)
db$neoadjuvant_therapy.[naj_na] <- db$neoadjuvant_therapy._2[naj_na]
db$mitotic_count_on_surgical

bio_idx <- db$Is_it_a_biopsy. == "yes"
sur_idx <- db$Is_it_a_biopsy. == "no"
 
str(db)

db_bio <- db[bio_idx,c(1,4,5,7,8)]
db_sur <- db[sur_idx,c(1,14,15)]
str(db_bio)
str(db_sur)

db <- merge(db_bio,db_sur)
colnames(db) <- c("PID", "Site", "Size", "m_bio", "Su", "resp_NA", "m_sur")
db <- db[complete.cases(db),]
db$Site <- factor(db$Site)
db$resp_NA <- factor(db$resp_NA)
db$Su <- as.numeric(db$Su)
summary(db[,-1])
pairs(db[,-1])
PD_idx <- db$m_bio <= db$m_sur & db$resp_NA == "yes" # identify the cases that m_sur > m_bio and have done therapy i.e. the non-responder
db$resp_NA[PD_idx] <- "no"
summary(db[,-1])
db$Site <- as.integer(db$Site)
db$resp_NA <- as.integer(db$resp_NA)
db <- db[db$resp_NA==1,-1]
db

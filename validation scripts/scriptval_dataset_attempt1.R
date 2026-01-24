# adjusting validation dataset for validation script

# libraries
library(dplyr)

# read the validation set file
validationdataset <- read.csv("validationset.csv")

# see all columns
colnames(validationdataset)

# keeping only interested columns 
datasetcol <- validationdataset%>% select(2, 5, 6, 8, 9, 15, 16)

# # check for duplicates in validationdataset
# duplicated(validationdataset)
# duplicated(validationdataset$Insert.Patient.ID)
# which(duplicated(validationdataset$Insert.Patient.ID))

# renaming columns to match simulation daset column names
colnames(datasetcol)[colnames(datasetcol) == "tumour.dimension.in.mm"] <- "size_mm"
colnames(datasetcol)[colnames(datasetcol) == "tumour.site"] <- "location"
colnames(datasetcol)[colnames(datasetcol) == "mitotic.count.on.biopsy"] <- "biopsy_mitosis"
colnames(datasetcol)[colnames(datasetcol) == "available.surface.on.biopsy..in.mm.2."] <- "biopsy_surface_hpf"
colnames(datasetcol)[colnames(datasetcol) == "mitotic.count.on.surgical"] <- "surgery_mitosis"
colnames(datasetcol)[colnames(datasetcol) == "neoadjuvant.therapy."] <- "neoadjuvant_therapy"

# change location names into numbers
datasetcol <- datasetcol %>%
  mutate(
    site = tolower(trimws(location)),
    location_code = case_when(
      grepl("stomach|gastric|gastro", site) ~ 4L,
      grepl("duoden", site) ~ 2L,
      grepl("small|ileum|jejunum|intestine", site) ~ 3L,
      grepl("colon|rectum|rectal|colorectal", site) ~ 1L,
      grepl("esophag|oesophag", site) ~ 5L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(-site)   # drop helper column

# drop location column
datasetcol$location <- NULL

# change column name
colnames(datasetcol)[colnames(datasetcol) == "location_code"] <- "location"

# convert biopsy surface from mm² to HPF
biopsy_surface_hpf <- as.numeric(
  datasetcol[["available.surface.on.biopsy..in.mm.2."]]
) * (23.5 / 5)

# change commas into dots
datasetcol$biopsy_surface_hpf <- as.numeric(gsub(",", ".", datasetcol$biopsy_surface_hpf))

# change from <1 to 0
datasetcol$biopsy_mitosis <- as.numeric(gsub("<1", "0", datasetcol$biopsy_mitosis))
datasetcol$surgery_mitosis <- as.numeric(gsub("<1", "0", datasetcol$surgery_mitosis))

# change yes vs no into 1 vs 0
datasetcol$neoadjuvant_therapy <- ifelse(datasetcol$neoadjuvant_therapy == "Yes", 1, 0)

# drop NA rows in size column since they're same patient id duplicate
datasetcol <- datasetcol[!is.na(datasetcol$size_mm), ]

# drop patient id column
datasetcol$Insert.Patient.ID <- NULL

# save modified dataset
write.csv(datasetcol, "validationsetmod.csv", row.names = FALSE)


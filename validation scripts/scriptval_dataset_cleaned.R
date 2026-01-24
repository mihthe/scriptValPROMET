# modifying validation dataset og

# load dataset
vset <- read.csv("validationsetog.csv")

# drop empty column
vset$Colonna.20 <- NULL

# make a working copy
vset_clean <- vset

# empty strings to NA
vset_clean[vset_clean == ""] <- NA

# from comma to period for surface columns
vset_clean$available.surface.on.biopsy..in.mm.2. <- gsub(
  ",", ".", vset_clean$available.surface.on.biopsy..in.mm.2.
)
vset_clean$available.surface.on.biopsy..in.mm.2..1 <- gsub(
  ",", ".", vset_clean$available.surface.on.biopsy..in.mm.2..1
)

# "<1" values to 0.5 (or 0?)
vset_clean$mitotic.count.on.biopsy <- gsub(
  "^<1$", "0.5", vset_clean$mitotic.count.on.biopsy
)
vset_clean$mitotic.count.on.surgical <- gsub(
  "^<1$", "0.5", vset_clean$mitotic.count.on.surgical
)


# ============================================
# compare mitotic count on biopsy vs 1 column
# ============================================

# filter biopsy rows
biopsy_rows <- vset_clean[vset_clean$Is.it.a.biopsy. == "yes", ]

# count non-NA values in both columns
main_non_na <- sum(!is.na(biopsy_rows$mitotic.count.on.biopsy))
dot1_non_na <- sum(!is.na(biopsy_rows$mitotic.count.on.biopsy.1))

# find rows with data in both columns, only in main column, and only in .1 column
both_have_data <- !is.na(biopsy_rows$mitotic.count.on.biopsy) & !is.na(biopsy_rows$mitotic.count.on.biopsy.1)
only_main <- !is.na(biopsy_rows$mitotic.count.on.biopsy) & is.na(biopsy_rows$mitotic.count.on.biopsy.1)
only_dot1 <- is.na(biopsy_rows$mitotic.count.on.biopsy) & !is.na(biopsy_rows$mitotic.count.on.biopsy.1)

# print results
message("\n=== BIOPSY MITOTIC COUNT ===\n")
message("Main column (mitotic.count.on.biopsy): Non-NA values: ", main_non_na)
message(".1 column (mitotic.count.on.biopsy.1): Non-NA values: ", dot1_non_na)
message("\nRows with data in BOTH columns: ", sum(both_have_data))
message("Rows with data ONLY in main column: ", sum(only_main))
message("Rows with data ONLY in .1 column: ", sum(only_dot1))

# compare values where both columns have data
if (sum(both_have_data) > 0) {
  comparison <- data.frame(
    main_col = biopsy_rows$mitotic.count.on.biopsy[both_have_data],
    dot1_col = biopsy_rows$mitotic.count.on.biopsy.1[both_have_data]
  )
  print(comparison)
}


# ============================================
# compare mitotic count on surgical vs 1 column
# ============================================

# filter surgical rows 
surgical_rows <- vset_clean[vset_clean$Is.it.a.biopsy. == "no", ]

# count non-NA values in both columns
main_non_na <- sum(!is.na(surgical_rows$mitotic.count.on.surgical))
dot1_non_na <- sum(!is.na(surgical_rows$mitotic.count.on.surgical.1))

# find rows with data in both columns, only in main column, and only in .1 column
both_have_surg <- !is.na(surgical_rows$mitotic.count.on.surgical) & !is.na(surgical_rows$mitotic.count.on.surgical.1)
only_main_surg <- !is.na(surgical_rows$mitotic.count.on.surgical) & is.na(surgical_rows$mitotic.count.on.surgical.1)
only_dot1_surg <- is.na(surgical_rows$mitotic.count.on.surgical) & !is.na(surgical_rows$mitotic.count.on.surgical.1)

# print results
message("\n=== SURGICAL MITOTIC COUNT ===\n")
message("Main column (mitotic.count.on.surgical): Non-NA values:", main_non_na)
message(".1 column (mitotic.count.on.surgical.1): Non-NA values:", dot1_non_na)
message("\nRows with data in BOTH columns:", sum(both_have_surg))
message("Rows with data ONLY in main column:", sum(only_main_surg))
message("Rows with data ONLY in .1 column:", sum(only_dot1_surg))

# compare values where both columns have data
if (sum(both_have_surg) > 0) {
  comparison_surg <- data.frame(
    main_col = surgical_rows$mitotic.count.on.surgical[both_have_surg],
    dot1_col = surgical_rows$mitotic.count.on.surgical.1[both_have_surg]
  )
  print(comparison_surg)
}


# ============================================
# compare available surface vs 1 column
# ============================================

# count non-NA values in both columns
main_non_na <- sum(!is.na(biopsy_rows$available.surface.on.biopsy..in.mm.2.))
dot1_non_na <- sum(!is.na(biopsy_rows$available.surface.on.biopsy..in.mm.2..1))

# find rows with data in both columns, only in main column, and only in .1 column
both_have_surf <- !is.na(biopsy_rows$available.surface.on.biopsy..in.mm.2.) & 
  !is.na(biopsy_rows$available.surface.on.biopsy..in.mm.2..1)
only_main_surf <- !is.na(biopsy_rows$available.surface.on.biopsy..in.mm.2.) & 
  is.na(biopsy_rows$available.surface.on.biopsy..in.mm.2..1)
only_dot1_surf <- is.na(biopsy_rows$available.surface.on.biopsy..in.mm.2.) & 
  !is.na(biopsy_rows$available.surface.on.biopsy..in.mm.2..1)

# print results
message("\n=== BIOPSY SURFACE ===\n")
message("Main column (available.surface.on.biopsy..in.mm.2.): Non-NA values:", main_non_na)
message(".1 column (available.surface.on.biopsy..in.mm.2..1): Non-NA values:", dot1_non_na)
message("\nRows with data in BOTH columns:", sum(both_have_surf))
message("Rows with data ONLY in main column:", sum(only_main_surf))
message("Rows with data ONLY in .1 column:", sum(only_dot1_surf))

# ============================================
# compare neoadjuvant vs 1 column
# ============================================
message("\n=== NEOADJUVANT THERAPY ===\n")
message("Main column (neoadjuvant.therapy.):Non-NA values:", sum(!is.na(surgical_rows$neoadjuvant.therapy.)))
message("\n.1 column (neoadjuvant.therapy..1):Non-NA values:", sum(!is.na(surgical_rows$neoadjuvant.therapy..1)))

both_have_nac <- !is.na(surgical_rows$neoadjuvant.therapy.) & 
  !is.na(surgical_rows$neoadjuvant.therapy..1)
message("\nRows with data in BOTH columns:", sum(both_have_nac))

# if overlap exists, check if values match
if (sum(both_have_nac) > 0) {
  cat("\nComparing values where BOTH exist:\n")
  comparison_nac <- data.frame(
    main = surgical_rows$neoadjuvant.therapy.[both_have_nac],
    dot1 = surgical_rows$neoadjuvant.therapy..1[both_have_nac]
  )
  cat("Values match:", sum(comparison_nac$main == comparison_nac$dot1), "/", nrow(comparison_nac), "\n")
}


# ============================================
# merge column and 1column
# ============================================

# new dataset
vset_merged <- vset_clean

# merge biopsy mitotic count columns
vset_merged$m_bio_combined <- ifelse(
  !is.na(vset_merged$mitotic.count.on.biopsy),
  vset_merged$mitotic.count.on.biopsy,
  as.character(vset_merged$mitotic.count.on.biopsy.1)
)

# merge surgical mitotic count columns
# correct "11/23" value in .1 column
vset_merged$mitotic.count.on.surgical.1 <- gsub(
  "11/23", "11", vset_merged$mitotic.count.on.surgical.1
)

vset_merged$m_surg_combined <- ifelse(
  !is.na(vset_merged$mitotic.count.on.surgical),
  vset_merged$mitotic.count.on.surgical,
  vset_merged$mitotic.count.on.surgical.1
)

# merge biopsy surface columns
# from comma to period for .1 column (as done before)
vset_merged$available.surface.on.biopsy..in.mm.2..1 <- gsub(
  ",", ".", vset_merged$available.surface.on.biopsy..in.mm.2..1
)

vset_merged$surface_combined <- ifelse(
  !is.na(vset_merged$available.surface.on.biopsy..in.mm.2.),
  vset_merged$available.surface.on.biopsy..in.mm.2.,
  vset_merged$available.surface.on.biopsy..in.mm.2..1
)

# merge neoadjuvant columns 
vset_merged$neoadjuvant_combined <- ifelse(
  !is.na(vset_merged$neoadjuvant.therapy.),
  vset_merged$neoadjuvant.therapy.,
  vset_merged$neoadjuvant.therapy..1
)

# ============================================
# identify paired patients
# ============================================
# patient IDs from biopsy rows
biopsy_patients <- vset_merged$Insert.Patient.ID[vset_merged$Is.it.a.biopsy. == "yes"]

# patient IDs from surgical rows
surgical_patients <- vset_merged$Insert.Patient.ID[vset_merged$Is.it.a.biopsy. == "no"]

# patients in BOTH biopsy and surgical rows
patients_in_both <- intersect(biopsy_patients, surgical_patients)

message("\n=== PATIENT PAIRING ===\n")
message("Unique patients with biopsy: ", length(unique(biopsy_patients)))
message("Unique patients with surgical: ", length(unique(surgical_patients)))
message("Patients appearing in BOTH: ", length(patients_in_both))

# ============================================
# create validation dataset with merged columns
# ============================================

# extract biopsy data for paired patients
biopsy_data <- vset_merged[vset_merged$Is.it.a.biopsy. == "yes" & 
                             vset_merged$Insert.Patient.ID %in% patients_in_both, ]

biopsy_clean <- data.frame(
  patient_id = biopsy_data$Insert.Patient.ID,
  site = biopsy_data$tumour.site,
  size_mm = biopsy_data$tumour.dimension.in.mm,
  m_bio = biopsy_data$m_bio_combined,
  surface_mm2 = biopsy_data$surface_combined
)

# extract surgical data for paired patients
surgical_data <- vset_merged[vset_merged$Is.it.a.biopsy. == "no" & 
                               vset_merged$Insert.Patient.ID %in% patients_in_both, ]

surgical_clean <- data.frame(
  patient_id = surgical_data$Insert.Patient.ID,
  m_surg = surgical_data$m_surg_combined,
  neoadjuvant = surgical_data$neoadjuvant_combined
)

# create dataset with merged columns
dat_validation <- merge(biopsy_clean, surgical_clean, by = "patient_id")

# convert to numeric
dat_validation$m_bio <- as.numeric(dat_validation$m_bio)
dat_validation$m_surg <- as.numeric(dat_validation$m_surg)
dat_validation$size_mm <- as.numeric(dat_validation$size_mm)
dat_validation$surface_mm2 <- as.numeric(dat_validation$surface_mm2)


# ============================================
# create validation dataset with "no" neoadjuvant
# ============================================

# option 1: exclude only confirmed "yes"
# option 2: keep only confirmed "no" (most conservative)
dat_option2 <- dat_validation[
  !is.na(dat_validation$neoadjuvant) & dat_validation$neoadjuvant == "no", 
]
cat("Option 2 (keep only 'no'):", nrow(dat_option2), "patients\n")

# ============================================
# create final validation dataset with complete cases
# ============================================

# find cases with ALL 5 variables present
complete <- complete.cases(
  dat_option2[, c("site", "size_mm", "m_bio", "surface_mm2", "m_surg")]
)
message("COMPLETE cases for validation:", sum(complete))

# create dataset with only complete cases
dat_final <- dat_option2[complete, ]
message("Final sample size:", nrow(dat_final))

# ============================================
# summary statistics of final validation dataset with complete cases
# ============================================

# numeric variables
message("Numeric variables:")
summary(dat_final[, c("size_mm", "m_bio", "surface_mm2", "m_surg")])

# site distribution
message("Tumor site distribution:")
print(table(dat_final$site))

# mean and median
cat("\nTumor size (mm):\n")
cat("  Mean (SD):", round(mean(dat_final$size_mm), 1), 
    "(", round(sd(dat_final$size_mm), 1), ")\n")
cat("  Median [IQR]:", median(dat_final$size_mm), 
    "[", quantile(dat_final$size_mm, 0.25), "-", 
    quantile(dat_final$size_mm, 0.75), "]\n")

cat("\nBiopsy mitotic count:\n")
cat("  Mean (SD):", round(mean(dat_final$m_bio), 2), 
    "(", round(sd(dat_final$m_bio), 2), ")\n")
cat("  Median [IQR]:", median(dat_final$m_bio), 
    "[", quantile(dat_final$m_bio, 0.25), "-", 
    quantile(dat_final$m_bio, 0.75), "]\n")

cat("\nSurgical mitotic count:\n")
cat("  Mean (SD):", round(mean(dat_final$m_surg), 2), 
    "(", round(sd(dat_final$m_surg), 2), ")\n")
cat("  Median [IQR]:", median(dat_final$m_surg), 
    "[", quantile(dat_final$m_surg, 0.25), "-", 
    quantile(dat_final$m_surg, 0.75), "]\n")

cat("\nBiopsy surface (mm²):\n")
cat("  Mean (SD):", round(mean(dat_final$surface_mm2), 2), 
    "(", round(sd(dat_final$surface_mm2), 2), ")\n")
cat("  Median [IQR]:", median(dat_final$surface_mm2), 
    "[", quantile(dat_final$surface_mm2, 0.25), "-", 
    quantile(dat_final$surface_mm2, 0.75), "]\n")

# ============================================
# location encoding
# ============================================

# tumours sites upon the Prometheus paper and app code:
# L = 1: Colonrectum
# L = 2: Duodenum  
# L = 3: Small intestine
# L = 4: Stomach

# create the Location variable
dat_final$L <- NA
dat_final$L[dat_final$site == "colonrectum"] <- 1
dat_final$L[dat_final$site == "duodenum"] <- 2
dat_final$L[dat_final$site == "small intestine"] <- 3
dat_final$L[dat_final$site == "stomach"] <- 4

# check the mapping
print(paste("Site to L mapping:"))
print(table(dat_final$site, dat_final$L, useNA = "ifany"))

# verify no missing L values
print(paste("\nMissing L values:", sum(is.na(dat_final$L)), "\n"))

# ============================================
# save final dataset
# ============================================
write.csv(dat_final, "validation_dataset_final.csv", row.names = FALSE)
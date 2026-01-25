# validation dataset from fulvia

# library
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(openxlsx)

# load dataset
valset <- read_excel("Validation PROMETEUS Renne.xlsx")

# drop empty columns
valset <- valset %>% select(-"Colonna 20")
valset <- valset %>% select(-"mitotic count on surgical 2")
valset <- valset %>% select(-"available surface on biopsy (in mm^2) 2")
valset <- valset %>% select(-"mitotic count on biopsy  2")
valset <- valset %>% select(-"digital?")
valset <- valset %>% select(-"digital? 2")

# drop column
valset <- valset %>% select(-"Informazioni cronologiche")

# drop therapy column
valset <- valset %>% select(-"neoadjuvant therapy?")
valset <- valset %>% select(-"neoadjuvant therapy? 2")

# remove rows due to missing data
valset <- valset[-c(5, 6, 11, 12, 35, 36, 41, 42, 45, 46, 87, 88, 91, 92, 93, 94, 99, 100, 103, 104), ]

# replace all "<1" with "0.5" (still text)
valset$`mitotic count on biopsy`[valset$`mitotic count on biopsy` == "<1"] <- "0.5"
valset$`mitotic count on surgical`[valset$`mitotic count on surgical` == "<1"] <- "0.5"

# change all "old" to "2821" in column '28.21'
valset$`available surface on biopsy (in mm^2)`[valset$`available surface on biopsy (in mm^2)` == "2821"] <- "28.21"

# save dataset
#write_xlsx(valset, "validation_set_lesscols.xlsx")
write_xlsx(valset, "/Users/mt/Documents/R_main/validation_set_lesscols.xlsx")

##


prepare_validation_data <- function(filepath, sheet = 1) {
  
  # read Excel safely
  raw <- read_excel(filepath, sheet = sheet, col_types = "text")
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  
  # expected columns (by position) 
  expected_cols <- c(
    "patient_id", "is_biopsy", "biopsy_N", "tumour_site", 
    "tumour_size_mm", "seen_SLR", "mitotic_biopsy", 
    "surface_biopsy_mm2", "surgical_N",
    "surgical_seen_SLR", "mitotic_surgery"
  )
  
  if (ncol(raw) < length(expected_cols)) {
    stop(
      sprintf(
        "Expected at least %d columns, but found %d.\nCheck the Excel file structure.",
        length(expected_cols), ncol(raw)
      )
    )
  }
  
  # rename only first expected columns (extra columns are ignored)
  colnames(raw)[1:length(expected_cols)] <- expected_cols
  
  # separate biopsy and surgery rows
  biopsy_rows <- raw %>%
    filter(tolower(is_biopsy) == "yes") %>%
    select(
      patient_id,
      tumour_site,
      tumour_size_mm,
      m_bio = mitotic_biopsy,
      surface_mm2 = surface_biopsy_mm2
    )
  
  surgery_rows <- raw %>%
    filter(tolower(is_biopsy) == "no") %>%
    select(
      patient_id,
      m_surg = mitotic_surgery,
    )
  
  # merge biopsy and surgery
  merged <- inner_join(biopsy_rows, surgery_rows, by = "patient_id")
  
  # data cleaning and transformations 
  merged <- merged %>%
    mutate(
      
      # mitotic count on biopsy
      m_bio = case_when(
        is.na(m_bio) ~ NA_real_,
        m_bio == "<1" ~ 0,
        TRUE ~ suppressWarnings(as.numeric(m_bio))
      ),
      
      # mitotic count on surgery
      m_surg = suppressWarnings(as.numeric(m_surg)),
      
      # tumour size
      size_mm = suppressWarnings(as.numeric(tumour_size_mm)),
      
      # surface area
      surface_mm2 = suppressWarnings(as.numeric(surface_mm2)),
      biopsy_surface_hpf = surface_mm2 / 0.2128,
      
      # tumour location mapping
      L = case_when(
        tolower(tumour_site) %in% c("colonrectum", "colon-rectum") ~ 1L,
        tolower(tumour_site) == "duodenum" ~ 2L,
        tolower(tumour_site) %in% c("small intestine", "small_intestine", "ileum") ~ 3L,
        tolower(tumour_site) %in% c("stomach", "gastric") ~ 4L,
        TRUE ~ NA_integer_
      ),
      
      # PROMETheus validation defaults
      n = 1L,
      y = 0L
    )
  
  # filter complete cases
  validation_ready <- merged %>%
    filter(
      !is.na(m_bio),
      !is.na(m_surg),
      !is.na(size_mm),
      !is.na(L),
      !is.na(surface_mm2)
    ) %>%
    select(
      patient_id,
      size_mm,
      L,
      m_bio,
      surface_mm2,
      m_surg,
      biopsy_surface_hpf,
      n,
      y
    )
  
  # reporting 
  cat(sprintf("Original rows in Excel: %d\n", nrow(raw)))
  cat(sprintf("Patients with biopsy + surgery: %d\n", nrow(merged)))
  cat(sprintf("Complete cases for validation: %d\n", nrow(validation_ready)))
  
  return(validation_ready)
}

# call the function using the dataset 
ready_data <- prepare_validation_data("validation_set_lesscols.xlsx")

# save the dataset
write.xlsx(ready_data, "validation_dataset_ready.xlsx")
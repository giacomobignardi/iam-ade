# SCRIPT:          tidy.R
# AUTHOR:          G.B. (Giaco)
# ADAPTED FROM:    Bignardi et al. (2022, Sci Rep), https://doi.org/10.1038/s41598-022-07161-z
# PURPOSE: prepare NTR data for iAM-ADE (DATE) modeling

rm(list = ls())

library(tidyverse) # Data wrangling
library(haven)     # Read .sav

# WORKING DIRECTORIES ####
# set open access working directories
wd_oa <- getwd()           # Local
wd_oa_scripts <- "R"       # Subdirectory for R scripts / functions
wd_oa_figures <- "Figures" # Subdirectory for output figures

# Secondary (non-open-access) working directory 
wd_noa <- substr(
  getwd(),
  0,
  nchar(getwd()) - nchar("analysis") - 1
)
wd_noa_data <- "working_data" # Subdirectory for data files

# DATA ####
# load twin data obtained from the NTR
dat <- read_sav(sprintf("%s/%s/upload_20260109_LL/PHE_20260109_5076_LL.sav", wd_noa, wd_noa_data))

# DESCRPITIVES ####
# variable names
colnames(dat)

# n indiviudals (per sex)
table(dat$sex)

# range age
range(dat$age7, na.rm = T)
range(dat$age8, na.rm = T)
range(dat$age10, na.rm = T)

# n families
dat %>%
  distinct(FID) %>%
  nrow()

# DATA SELECTION ####
# select family members
dat_estd <- dat %>%
  dplyr::filter(Extension == 1 | Extension == 2 | Extension == 10 | Extension == 61 | Extension == 62) %>% # select twins
  # 1 first multiple in family
  # 2 2nd multiple in family
  # 10 sibling
  # 61 first registered spouse of person 1 in family
  # 62 first registered spouse of person 2 in family
  dplyr::filter(multiple_type <= 2) %>% # filter for < triplet
  dplyr::select(c(FID, Extension, twzyg, sex, 
                  age7, age8, age10, 
                  starts_with(c("NEO", "neo", "BOR", "bor", "len","prtjr", "ea")),
                  )) %>%
  dplyr::mutate(across(.cols = starts_with(c("NEO", "neo", "BOR", "bor", "len","prtjr", "ea")),
                       .fns = as.numeric)) %>%
  #   value label
  # 1   MZM
  # 2   DZM
  # 3   MZF
  # 4   DZF
  # 5 DOSmf
  # 6 DOSfm
  dplyr::filter(twzyg == 1 | twzyg == 2 | twzyg == 3 | twzyg == 4 | twzyg == 5 | twzyg == 6 | is.na(twzyg)) %>%
  # 1 male 2 female
  dplyr::filter(sex == 1 | sex == 2)


# SAVE --------------------------------------------------------------
#family twin data ####
write.csv(as.data.frame(dat_estd),sprintf("%s/%s/01_dat_estd.csv",wd_noa,wd_noa_data))

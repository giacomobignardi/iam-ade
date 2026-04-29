# SCRIPT:  tidy.R
# AUTHOR:  G.B. (Giaco)
# PURPOSE: Estimate between-individual variance (ICC) using repeated measures

rm(list = ls())

library(tidyverse) # Data wrangling
library(psych)     # for ICC

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

# REPRODUCIBILITY ####
set.seed(42)

# DATA ####
# load selected fam of twin data
dat <- read.csv(sprintf("%s/%s/01_dat_estd.csv", wd_noa, wd_noa_data))
set.seed(42)

# remove individuals of age < 18
dat <- dat %>%
  filter((age7 >= 18 | is.na(age7)) &
           (age8 >= 18 | is.na(age8)) &
           (age10 >= 18 | is.na(age10)))

# make colnames lower case
colnames(dat) <- tolower(colnames(dat))

# rename bor sum score
dat <- dat %>% rename(bpd_7 = bor_7, bpd_8 = bor_8, bpd_10 = bor_10)
var_name <- "bpd"

# ICC  ####
rYY <- psych::ICC(dat %>% dplyr::select(starts_with(var_name)))$results %>% filter(type == "ICC2")
save(rYY, file = sprintf("%s/%s/04.01_rYY_bpd.Rdata", wd_noa, wd_noa_data))


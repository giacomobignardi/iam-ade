# SCRIPT:        sampling.R
# AUTHOR:        G.B. (Giaco)
# ADAPTED FROM:  Bignardi et al. (2022, Sci Rep), https://doi.org/10.1038/s41598-022-07161-z
# PURPOSE:       Pseudo-randomized (optimized) sampling (with twins as unit of sampling)

rm(list = ls())

library(tidyverse) # Data wrangling

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
# load selected fam of twin data
dat <- read.csv(sprintf("%s/%s/01_dat_estd.csv", wd_noa, wd_noa_data))
set.seed(42)
# whole sample descriptive
nrow(dat %>% distinct(FID))
# 10296 families
nrow(dat)
# 19346 total individuals
nrow(dat %>% filter(sex == 2))
# 12315 women
dat %>% reframe(range(age7, na.rm = T))
dat %>% reframe(range(age8, na.rm = T))
dat %>% reframe(range(age10, na.rm = T))
# > dat %>% reframe(range(age8, na.rm = T))
# range(age8, na.rm = T)
# 1                     11
# 2                     97

# make colnames lower case
tolower(colnames(dat))
colnames(dat) <- tolower(colnames(dat))

# add space between var name and var wave
var_name <- "ea7_agg"

# 1 first multiple in family
# 2 2nd multiple in family
# 10 sibling
# 61 first registered spouse of person 1 in family
# 62 first registered spouse of person 2 in family

# SAMPLING ####
# sampling procedure
# select twins and make a wide dataframe with one row per twin pair across waves of data collection
dat_tw <- dat %>%
  filter(!is.na(twzyg)) %>% # keep only twins with known zygosity
  rename(age = ea7_age_agg) %>%
  dplyr::select(c(fid, extension, twzyg, age, var_name)) %>%
  pivot_wider(names_from = "extension", values_from = "age":var_name) %>%
  rename(
    age_tw1 = age_1, age_tw2 = age_2,
    tw1 = paste0(var_name,"_1"), tw2 = paste0(var_name,"_2")
  ) %>%
  mutate( # add sex (to use as a covariate later MALE 1, FEMALE 2)
    sex_tw1 = ifelse(twzyg == 1 | twzyg == 2 | twzyg == 5, 1, 2),
    sex_tw2 = ifelse(twzyg == 1 | twzyg == 2 | twzyg == 6, 1, 2)
  )
# final sample
table(dat_tw$twzyg)
# 1    2    3    4    5    6 
# 1303  887 2867 1758 1284 1275 

# simplify descriptives as we do not need information of opposite sex twins
dat_tw$twzyg <- ifelse(dat_tw$twzyg == 1 | dat_tw$twzyg == 3, "MZ", "DZ")

# descriptive twins
tw_des <- dat_tw %>%
  group_by(twzyg) %>%
  reframe(
    n = sum(!is.na(tw1)),
    mean = mean(age_tw1, na.rm = T),
    sd = sd(age_tw1, na.rm = T)
  )

# extend to family (sibs,spouses)
dat_fam <-
  dat %>%
  dplyr::filter(is.na(twzyg) & !(extension == 1 | extension == 2)) %>%
  rename(age = ea7_age_agg) %>%
  dplyr::select(c(fid, sex, extension, age, var_name)) %>%
  pivot_wider(names_from = "extension", values_from = c("sex", "age":var_name)) %>%
  rename(
    age_s1 = age_10, age_stw1 = age_61, age_stw2 = age_62,
    sex_s1 = sex_10, sex_stw1 = sex_61, sex_stw2 = sex_62,
    s1 = paste0(var_name,"_10"), stw1 = paste0(var_name,"_61"), stw2 = paste0(var_name,"_62")
  ) %>%
  semi_join(dat_tw, by = c("fid"))  # match to wave sampled for twins (twins are sampling units)

sum(!is.na(dat_fam$s1))
#[1] 1233
sum(!is.na(dat_fam$stw1))
#[1] 660
sum(!is.na(dat_fam$stw2))
#[1] 651

# create final extended twin family including partners
dat_twfam <- merge(dat_tw, dat_fam, c("fid"), all.x = T) %>%
  arrange(fid) %>%
  # årrange columns to be easier to read
  dplyr::select(
    fid,
    twzyg,
    age_tw1, age_tw2, age_s1, age_stw1, age_stw2,
    sex_tw1, sex_tw2, sex_s1, sex_stw1, sex_stw2,
    tw1, tw2, s1, stw1, stw2,
  ) 


# remove same sex partners
nrow(dat_twfam %>% filter(sex_tw1 == sex_stw1))
# 22
nrow(dat_twfam %>% filter(sex_tw2 == sex_stw2))
# 13
dat_twfam[which(dat_twfam$sex_tw1 == dat_twfam$sex_stw1), ]$stw1  <- NA
dat_twfam[which(dat_twfam$sex_tw2 == dat_twfam$sex_stw2), ]$stw2  <- NA

# SAMPLE ####
# total sample
sum(!is.na(dat_twfam$tw1)) +
  sum(!is.na(dat_twfam$tw2)) +
  sum(!is.na(dat_twfam$s1)) +
  sum(!is.na(dat_twfam$stw1)) +
  sum(!is.na(dat_twfam$stw2))
#[1] 11206

# n of women
table(dat_twfam[which(!is.na(dat_twfam$tw1)), ]$sex_tw1) +
  table(dat_twfam[which(!is.na(dat_twfam$tw2)), ]$sex_tw2) +
  table(dat_twfam[which(!is.na(dat_twfam$s1)), ]$sex_s1) +
  table(dat_twfam[which(!is.na(dat_twfam$stw1)), ]$sex_stw1) +
  table(dat_twfam[which(!is.na(dat_twfam$stw2)), ]$sex_stw2)
# 1     2 
# 3968 7238

# n of families
nrow(dat_twfam %>% distinct(fid))
#[1] 9374

# check that there is no individual younger then 18
range(dat_twfam$age_tw1, na.rm = T)
range(dat_twfam$age_tw2, na.rm = T)
range(dat_twfam$age_s1, na.rm = T)
range(dat_twfam$age_stw1, na.rm = T)
range(dat_twfam$age_stw2, na.rm = T)

# MZ
n_MZ_m <- table(dat_twfam[which(!is.na(dat_twfam$tw1)), ][dat_twfam$twzyg == "MZ", ]$sex_tw1)[1] + table(dat_twfam[which(!is.na(dat_twfam$tw2)), ][dat_twfam$twzyg == "MZ", ]$sex_tw2)[1]
n_MZ_f <- table(dat_twfam[which(!is.na(dat_twfam$tw1)), ][dat_twfam$twzyg == "MZ", ]$sex_tw1)[2] + table(dat_twfam[which(!is.na(dat_twfam$tw2)), ][dat_twfam$twzyg == "MZ", ]$sex_tw2)[2]
# DZ
n_DZ_m <- table(dat_twfam[which(!is.na(dat_twfam$tw1)), ][dat_twfam$twzyg == "DZ", ]$sex_tw1)[1] + table(dat_twfam[which(!is.na(dat_twfam$tw2)), ][dat_twfam$twzyg == "DZ", ]$sex_tw2)[1]
n_DZ_f <- table(dat_twfam[which(!is.na(dat_twfam$tw1)), ][dat_twfam$twzyg == "DZ", ]$sex_tw1)[2] + table(dat_twfam[which(!is.na(dat_twfam$tw2)), ][dat_twfam$twzyg == "DZ", ]$sex_tw2)[2]
# Full-sibs
n_s_m <- table(dat_twfam[which(!is.na(dat_twfam$s1)), ]$sex_s1)[1]
n_s_f <- table(dat_twfam[which(!is.na(dat_twfam$s1)), ]$sex_s1)[2]
# husbands
n_h <- table(dat_twfam[which(!is.na(dat_twfam$stw1)), ]$sex_stw1)[1] + table(dat_twfam[which(!is.na(dat_twfam$stw2)), ]$sex_stw2)[1]
# wifes
n_w <- table(dat_twfam[which(!is.na(dat_twfam$stw1)), ]$sex_stw1)[2] + table(dat_twfam[which(!is.na(dat_twfam$stw2)), ]$sex_stw2)[2]

# final sample size
dem <- data.frame(
  N = c(n_MZ_m, n_MZ_f, n_DZ_m, n_DZ_f, n_s_m, n_s_f, n_h, n_w),
  Type = c("MZ", "MZ", "DZ", "DZ", "Full sib", "Full sib", "partner", "partner"),
  Sex = c("men", "women", "men", "women", "men", "women", "man", "women")
)

# SAVE ####
# demographic
write_csv(dem, sprintf("%s/SI/SI.01_table1_dem_EA.csv", wd_oa))

## sampled family twin data ####
write_csv(dat_twfam, sprintf("%s/%s/SI.01_dat_estd_EA.csv", wd_noa, wd_noa_data))


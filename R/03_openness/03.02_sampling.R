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
dat <- dat %>% rename(oe_7 = neo_o7, oe_8 = neo_o8, oe_10 = neo_o10)
var_name <- "oe"

# sample descriptives for individuals with at least one measure across surveys
dat_complete <- dat %>% filter(!(is.na(paste0(var_name,"_7")) & is.na(paste0(var_name,"_8")) & is.na(paste0(var_name,"_10"))))
nrow(dat_complete)
# 19346
nrow(dat_complete %>% filter(sex == 2))
# 12315

table(dat_complete$extension)
# 1    2   10   61   62 
# 7978 7529 2278  791  770

# 1 first multiple in family
# 2 2nd multiple in family
# 10 sibling
# 61 first registered spouse of person 1 in family
# 62 first registered spouse of person 2 in family

# SAMPLING --------------------------------------------------------------
# sampling procedure
# select twins and make a wide dataframe with one row per twin pair across waves of data collection
dat_tw <- dat %>%
  filter(!is.na(twzyg)) %>% # keep only twins with known zygosity
  rename(age_7 = age7, age_8 = age8, age_10 = age10) %>%
  dplyr::select(c(fid, extension, twzyg, age_7, age_8, age_10, paste0(var_name,"_7"), paste0(var_name,"_8"), paste0(var_name,"_10"))) %>%
  pivot_longer(age_7:(paste0(var_name,"_10")), names_to = c("var", "wave"), names_sep = "_", values_to = "score") %>%
  mutate(score = as.numeric(score)) %>%
  pivot_wider(names_from = c("var"), values_from = "score") %>%
  pivot_wider(names_from = "extension", values_from = "age":var_name) %>%
  rename(
    age_tw1 = age_1, age_tw2 = age_2,
    tw1 = paste0(var_name,"_1"), tw2 = paste0(var_name,"_2")
  ) %>%
  mutate(complete_data = is.na(tw1) + is.na(tw2)) # index for complete pair: 2 = both missing, 1 = one complete, 0 = both complete

# select pairs with complete data
comp_tw <- dat_tw %>%
  dplyr::filter(complete_data == 0) %>%
  group_by(fid) %>%
  slice_sample(n = 1) # sample only one per family across waves

# select pairs with data only for only one twin
sing_tw <- dat_tw %>%
  dplyr::filter(complete_data == 1) %>%
  group_by(fid) %>%
  slice_sample(n = 1) %>% # sample only one per family across waves
  filter(!(fid %in% comp_tw$fid))

# final data frame for twins
dat_twf <- rbind(comp_tw, sing_tw) %>%
  arrange(fid) %>%
  mutate( # add sex (to use as a covariate later MALE 1, FEMALE 2)
    sex_tw1 = ifelse(twzyg == 1 | twzyg == 2 | twzyg == 5, 1, 2),
    sex_tw2 = ifelse(twzyg == 1 | twzyg == 2 | twzyg == 6, 1, 2)
  ) %>%
  dplyr::select(-complete_data)

# final sample
table(dat_twf$twzyg)
# 1    2    3    4    5    6 
# 1254  849 2779 1709 1228 1214 

# impute missing age
# how many first-born twins have missing age? 0
sum(!is.na(dat_twf$tw1) & is.na(dat_twf$age_tw1))

# how many second-born twins have missing age? 0
sum(!is.na(dat_twf$tw2) & is.na(dat_twf$age_tw2))

# remove individual with no age information
nrow(dat_twf)
dat_twf <- dat_twf %>%
  filter(!(!is.na(tw1) & is.na(age_tw1))) %>%
  filter(!(!is.na(tw2) & is.na(age_tw2)))

nrow(dat_twf)
table(dat_twf$twzyg)
# there is no individual with missing age

# make age equal across twins
dat_twf <- dat_twf %>%
  mutate(
    age_tw1 = ifelse(is.na(age_tw1), age_tw2, age_tw1),
    age_tw2 = ifelse(is.na(age_tw2), age_tw1, age_tw2)
  ) 

# which twins report different ages
table(dat_twf$age_tw1-dat_twf$age_tw2) 
# -5   -4   -3   -2   -1    0    1    2    3    4 
# 2    4   58   66  283 8224  283   55   54    4  
nrow(dat_twf[which((dat_twf$age_tw1 - dat_twf$age_tw2) != 0),])
# for 809 twins age was different, with a maximum difference of 5 years apart

# remove individuals of age less then 18
nrow(dat_twf %>% filter(age_tw1 < 18 | age_tw2 < 18 ))
# 88
table(dat_twf[dat_twf$age_tw1 < 18 | dat_twf$age_tw2 < 18, ]$twzyg)
# 1  2  3  4  5  6 
# 15 12 33 17  7  4 

dat_twf <- dat_twf %>% filter(age_tw1 > 17 & age_tw2 > 17)
table(dat_twf$twzyg)
# 1    2    3    4    5    6 
# 1239  837 2746 1692 1221 1210 

# simplify descriptives as we do not need information of opposite sex twins
dat_twf$twzyg <- ifelse(dat_twf$twzyg == 1 | dat_twf$twzyg == 3, "MZ", "DZ")

# descriptive twins
tw_des <- dat_twf %>%
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
  rename(age_7 = age7, age_8 = age8, age_10 = age10,
         prtjr_8 = prtjr8 , prtjr_10 = prtjr10) %>%
  dplyr::select(c(fid, sex, extension, age_7, age_8, age_10, prtjr_8, prtjr_10, paste0(var_name,"_7"), paste0(var_name,"_8"), paste0(var_name,"_10"))) %>%
  pivot_longer(age_7:paste0(var_name,"_10"), names_to = c("var", "wave"), names_sep = "_", values_to = "score") %>%
  mutate(score = as.numeric(score)) %>%
  pivot_wider(names_from = c("var"), values_from = "score") %>%
  pivot_wider(names_from = "extension", values_from = c("sex", "age":var_name)) %>%
  rename(
    age_s1 = age_10, age_stw1 = age_61, age_stw2 = age_62,
    sex_s1 = sex_10, sex_stw1 = sex_61, sex_stw2 = sex_62,
    prtjr_stw1 = prtjr_61, prtjr_stw2 = prtjr_62,
    s1 = paste0(var_name,"_10"), stw1 = paste0(var_name,"_61"), stw2 = paste0(var_name,"_62")
  ) %>%
  semi_join(dat_twf, by = c("fid", "wave")) %>%  # match to wave sampled for twins (twins are sampling units)
  rowwise() %>%
  mutate(prtjr = ifelse(sum(!is.na(c(prtjr_stw1,prtjr_stw2))) ==1, 
                        coalesce(prtjr_stw1,prtjr_stw2),
                        rowMeans(cbind(prtjr_stw1,prtjr_stw2), na.rm = T))) %>% # obtain average years of relationship per pair
  dplyr::select(-c(prtjr_10, prtjr_stw1, prtjr_stw2))

# add family (sibs,spouses) with no twins
View(dat %>% dplyr::filter(!(extension == 1 | extension == 2) & !(fid %in% dat_twf$fid)) %>% 
  rename(age_7 = age7, age_8 = age8, age_10 = age10,
         prtjr_8 = prtjr8 , prtjr_10 = prtjr10) %>%
  dplyr::select(c(fid, sex, twzyg,extension, age_7, age_8, age_10, prtjr_8, prtjr_10, paste0(var_name,"_7"), paste0(var_name,"_8"), paste0(var_name,"_10"))) %>%
  pivot_longer(age_7:paste0(var_name,"_10"), names_to = c("var", "wave"), names_sep = "_", values_to = "score") %>%
  mutate(score = as.numeric(score)) %>%
  pivot_wider(names_from = c("var"), values_from = "score") %>%
  pivot_wider(names_from = "extension", values_from = c("sex", "age":var_name)) %>%
  rename(
    age_s1 = age_10, age_stw1 = age_61, age_stw2 = age_62,
    sex_s1 = sex_10, sex_stw1 = sex_61, sex_stw2 = sex_62,
    prtjr_stw1 = prtjr_61, prtjr_stw2 = prtjr_62,
    s1 = paste0(var_name,"_10"), stw1 = paste0(var_name,"_61"), stw2 = paste0(var_name,"_62")
  ))

sum(!is.na(dat_fam$s1))
#[1] 1330
sum(!is.na(dat_fam$stw1))
#[1] 504
sum(!is.na(dat_fam$stw2))
#[1] 486

# create final extended twin family including partners
dat_twfam <- merge(dat_twf, dat_fam, c("fid", "wave"), all.x = T) %>%
  arrange(fid) %>%
  # årrange columns to be easier to read
  dplyr::select(
    fid, wave,
    twzyg,
    age_tw1, age_tw2, age_s1, age_stw1, age_stw2,
    sex_tw1, sex_tw2, sex_s1, sex_stw1, sex_stw2,
    tw1, tw2, s1, stw1, stw2,
    prtjr
  ) 

# check for missing ages (here imputation is not possible)
sum(!is.na(dat_twfam$s1) & is.na(dat_twfam$age_s1))
# 0
sum(!is.na(dat_twfam$stw1) & is.na(dat_twfam$age_stw1))
# 0
sum(!is.na(dat_twfam$stw2) & is.na(dat_twfam$age_stw2))
# 0

# remove individuals of age less then 18
nrow(dat_twfam %>% filter(age_s1 < 18))
# 30
table(dat_twfam[dat_twfam$age_s1 < 18, ]$sex_s1)
# 1  2 
# 15 15 
dat_twfam[which(!is.na(dat_twfam$s1) & dat_twfam$age_s1 < 18), ]$s1 <- NA
dat_twfam[which(dat_twfam$age_s1 < 18), ]$age_s1 <- NA

# remove same sex partners
nrow(dat_twfam %>% filter(sex_tw1 == sex_stw1))
# 22
nrow(dat_twfam %>% filter(sex_tw2 == sex_stw2))
# 14
dat_twfam[which(dat_twfam$sex_tw1 == dat_twfam$sex_stw1), ]$stw1  <- NA
dat_twfam[which(dat_twfam$sex_tw2 == dat_twfam$sex_stw2), ]$stw2  <- NA

# SAMPLE --------------------------------------------------------------
# total sample
sum(!is.na(dat_twfam$tw1)) +
  sum(!is.na(dat_twfam$tw2)) +
  sum(!is.na(dat_twfam$s1)) +
  sum(!is.na(dat_twfam$stw1)) +
  sum(!is.na(dat_twfam$stw2))
#[1] 16235

# n of women
table(dat_twfam[which(!is.na(dat_twfam$tw1)), ]$sex_tw1) +
  table(dat_twfam[which(!is.na(dat_twfam$tw2)), ]$sex_tw2) +
  table(dat_twfam[which(!is.na(dat_twfam$s1)), ]$sex_s1) +
  table(dat_twfam[which(!is.na(dat_twfam$stw1)), ]$sex_stw1) +
  table(dat_twfam[which(!is.na(dat_twfam$stw2)), ]$sex_stw2)
# 1     2 
# 5710 10525

# n of families
nrow(dat_twfam %>% distinct(fid))
#[1] 8945

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

# SAVE --------------------------------------------------------------
# demographic
write_csv(dem, sprintf("%s/SI/03.02_table_dem_openness.csv", wd_oa))

# sampled family twin data ####
write_csv(dat_twfam, sprintf("%s/%s/03.02_dat_estd_openness.csv", wd_noa, wd_noa_data))

# SCRIPT:  proxy_conv.R
# AUTHOR:  G.B. (Giaco)
# PURPOSE: Test for convergence as alternative explanation for observed partner correlations

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

# REPRODUCIBILITY ####
set.seed(42)

# DATA ####
# load extended twin family data
dat_fam2g <- read.csv(sprintf("%s/%s/03.02_dat_estd_openness.csv", wd_noa, wd_noa_data))

# match men and women (always in the same column)
col_names <- c("par1", "par2", "sex_par1", "sex_par2", "len_rel")
dat_par <- rbind(
  dat_fam2g %>% filter(sex_stw1 == 1 & !is.na(stw1) & !is.na(tw1) ) %>% dplyr::select(stw1, tw1, sex_stw1,  sex_tw1, prtjr) %>% rename_with(~col_names),
  dat_fam2g %>% filter(sex_stw1 == 2 & !is.na(stw1) & !is.na(tw1) ) %>% dplyr::select(tw1, stw1, sex_tw1,  sex_stw1, prtjr) %>% rename_with(~col_names),
  dat_fam2g %>% filter(sex_stw2 == 1 & !is.na(stw2) & !is.na(tw2) ) %>% dplyr::select(stw2, tw2, sex_stw2,  sex_tw2, prtjr) %>% rename_with(~col_names),
  dat_fam2g %>% filter(sex_stw2 == 2 & !is.na(stw2) & !is.na(tw2) ) %>% dplyr::select(tw2, stw2, sex_tw2,  sex_stw2, prtjr) %>% rename_with(~col_names)
)

# does the relationship between partners depend upon the length of their relationshp?
sumy <- summary(lm(scale(par1)~scale(par2)*scale(len_rel), dat_par))
sumy <- sumy$coefficients %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%  
  filter(rowname == "scale(par2):scale(len_rel)") %>% 
  mutate(
    across(c("Estimate", "Std. Error","t value", "Pr(>|t|)"), ~ round(.x, 3))
  ) %>% 
  rename(
    estimate  = Estimate ,
    SE  = "Std. Error",
    t  = "t value",
    p = "Pr(>|t|)"
  ) %>% 
  select(-rowname)

#SAVE####
write_csv(sumy, sprintf("%s/SI/03.05_table_lm_conv_openness.csv", wd_oa))

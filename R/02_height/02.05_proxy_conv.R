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
dat_fam2g <- read.csv(sprintf("%s/%s/02.02_dat_estd_height.csv", wd_noa, wd_noa_data))

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
write_csv(sumy, sprintf("%s/SI/02.05_table_lm_conv_height.csv", wd_oa))


# Residuals:
#     Min      1Q  Median      3Q     Max 
# -2.8165 -0.5724 -0.0100  0.5573  2.5379 
# 
# Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.06358    0.04367   1.456  0.14614    
# scale(par2)                 0.16805    0.04396   3.823  0.00015 ***
# scale(len_rel)             -0.11054    0.04457  -2.480  0.01348 *  
# scale(par2):scale(len_rel)  0.02132    0.04349   0.490  0.62417    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9164 on 461 degrees of freedom
#   (595 observations deleted due to missingness)
# Multiple R-squared:  0.05626,	Adjusted R-squared:  0.05012 
# F-statistic: 9.161 on 3 and 461 DF,  p-value: 6.729e-06
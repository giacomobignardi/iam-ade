# SCRIPT:  iAMDATE.R
# AUTHOR:  G.B. (Giaco)
# PURPOSE: Apply iAM-ADE (DATE) model to estimate variance components accounting for AM.

rm(list = ls())

library(tidyverse)  # Data wrangling and plotting
library(OpenMx)     # Structural equation modelling

# WORKING DIRECTORIES ####
# Flag: set TRUE to use HPC cluster paths, FALSE for local development
is_cluster <- T

# Primary (open-access) working directory
wd_oa <- ifelse(
  is_cluster,
  "/data/workspaces/lag/workspaces/lg-iamdate/analysis",   # HPC (Linux/NPSOL)
  getwd()                                                  # Local
)

wd_oa_scripts <- "R"        # Subdirectory for R scripts / functions
wd_oa_figures <- "Figures"  # Subdirectory for output figures

# Secondary (non-open-access) working directory 
wd_noa <- substr(wd_oa, 0, nchar(wd_oa) - nchar("analysis") - 1)
wd_noa_data <- "working_data"  # Subdirectory for data files

# OPTIMIZER & REPRODUCIBILITY ####
mxOption(NULL, "Default optimizer", "NPSOL")   # Use NPSOL for OpenMx fitting
set.seed(42)                                   # Seed for reproducibility

# DATA ####
# load extended twin fadmily data
dat_fam2g <- read.csv(sprintf("%s/%s/SI.01_dat_estd_EA.csv", wd_noa, wd_noa_data))
dat_fam2g <- dat_fam2g %>% 
  rename(value_tw1 = tw1, 
         value_tw2 = tw2, 
         value_s1 = s1, 
         value_stw1 = stw1, 
         value_stw2 = stw2)

# make data long
dat_fam2g_long <- dat_fam2g %>% 
  pivot_longer(cols = age_tw1:value_stw2, names_to = c(".value", "rel"), names_sep = "_") %>% 
  # recode EA in years of education (for consistency with previous lit. Gonggrijp et al., 2023 DOI: 10.3389/fgene.2023.1150697)
  mutate(value = case_when(
    value == 1 ~ 6,
    value %in% c(2, 3) ~ 10,
    value %in% c(4, 5) ~ 13,
    value == 6 ~ 15,
    value == 7 ~ 16,
    TRUE ~ NA_integer_  # Optional: handle unexpected values
  )) %>% 
  # scale within sex (for consistency Sunde et al. 2025 https://doi.org/10.1038/s41467-025-60483-0)
  group_by(sex) %>% 
  mutate(z = as.numeric(scale(value))) %>% 
  ungroup() 

# FUNCTIONS & DATA PREPARATION ####
dat_fam2g <- dat_fam2g_long %>%
  dplyr::select(fid, twzyg, rel, z) %>%
  pivot_wider(names_from = rel, values_from = z) %>%
  rename(
    value_tw1  = tw1,
    value_tw2  = tw2,
    value_s1   = s1,
    value_stw1 = stw1,
    value_stw2 = stw2
  )

# Load models & source functions ####
# source ctd and estd inspired model
source(sprintf("%s/%s/function/iam.date.fun.R", wd_oa, wd_oa_scripts))
source(sprintf("%s/%s/function/acde.fun.R", wd_oa, wd_oa_scripts))

# select family members for analysis and adapt to specifics of iAM model
# remember that the order needs to match exactly the order in the function
vars <- "value"
sel_vars <- paste0("value", "_", c("stw1", "tw1", "s1", "tw2", "stw2"))
# select data for analysis
mz <- subset(dat_fam2g, twzyg == "MZ", sel_vars)
dz <- subset(dat_fam2g, twzyg == "DZ", sel_vars)
mz_sample <- psych::pairwiseCount(mz)
dz_sample <- psych::pairwiseCount(dz)

# name for readability
colnames(mz_sample) <- sel_vars
rownames(mz_sample) <- sel_vars
colnames(dz_sample) <- sel_vars
rownames(dz_sample) <- sel_vars
# source helper to fit estd inspierd model

# before loading source it is important to specify expected means
mean_MZ <- c("mean", "mean", "mean", "mean", "mean")
mean_DZ <- c("mean", "mean", "mean", "mean", "mean")

# ACE ####
# first fit classical twin design
# with common environement
ace_mod <- acde.fun(mz, dz, c("value_tw1", "value_tw2"), dp.free = F, cp.free = T)
ace_fit <- mxTryHard(ace_mod, intervals = T)
ace_fit$US$result

# without shared environmental effects
ae_mod <- acde.fun(mz, dz, c("value_tw1", "value_tw2"))
ae_fit <- mxTryHard(ae_mod, intervals = T)
ae_fit$US$result

# final models
ace_mxcompare <- as.data.frame(mxCompare(ace_fit, ae_fit)) %>%
  dplyr::select(-c(fitUnits, fit, diffFit, chisq, SBchisq)) %>%
  mutate(minus2LL = round(minus2LL, 0), AIC = round(AIC, 0), diffLL = round(diffLL, 2), p = round(p, 5))

ace_est <- summary(ace_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c("V","A","D","C","E"), model = "ACE")

# Save table CTD
write_csv(ace_est %>% select(model, parameter, estimate, lbound, ubound) %>% mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5)), sprintf("%s/SI/SI.04_tableCTD_est_EA.csv", wd_oa))

# HELPER FUNCTION ####
# Note that sourcing this will overwrite some model specification (e.g. ACE/ADE)
source(sprintf("%s/%s/function/iam.date.fun.helper.R", wd_oa, wd_oa_scripts))

# iAM-DATE --------------------------------------------------------------
iam_mods <- c(
  xiam_date_list, qiam_date_list, iam_date_list, date_list,
  xiam_cate_list, qiam_cate_list, iam_cate_list, cate_list
)

# list of mxModels to fit to data
iam_results <- lapply(iam_mods, run.model.safe)

# Extract diagnostic and AIC 
iam_results_diagnostic <- map_dfr(iam_results, ~ tibble(
  model_name   = .x$model_name,
  status_code  = .x$status_code,
  status_msg   = .x$status_msg,
  used_tryhard = .x$used_tryhard,
  AIC          = AIC(.x$fit),
)) %>%
  mutate(index = row_number()) %>%  
  mutate(delta_AIC = AIC - min(AIC)) 
# Extract diagnostic and AIC
iam_results_diagnostic %>% filter(status_code != 0)

# MODEL COMPARISON ####
# Indirect-direct assortment
# indirect assortment
# iam_cate_fit <- iam_results[[iam_results_diagnostic %>% filter(model_name == "iam_cate" ) %>% pull(index)]]$fit
# iam_ace_fit <- iam_results[[iam_results_diagnostic %>% filter(model_name == "iam_ace" ) %>% pull(index)]]$fit
# mxCompare(iam_cate_fit, iam_ace_fit)
# 
# dam_ace_fit <- iam_results[[iam_results_diagnostic %>% filter(model_name == "dam_ace" )%>% slice(1) %>% pull(index)]]$fit
# mxCompare(iam_ace_fit, dam_ace_fit)
# 
# iam_ace_fit$est_var 
# iam_ace_fit$est_path
# biologically plausible
qiam_cate_fit <- iam_results[[iam_results_diagnostic %>% filter(model_name == "qiam_cate" ) %>% pull(index)]]$fit
qiam_ace_fit  <- iam_results[[iam_results_diagnostic %>% filter(model_name == "qiam_ace" ) %>% pull(index)]]$fit
qdam_ace_fit  <- iam_results[[iam_results_diagnostic %>% filter(model_name == "qdam_ace" ) %>% pull(index)]]$fit
iam_ace_fit   <- iam_results[[iam_results_diagnostic %>% filter(model_name == "iam_ace" ) %>%slice(1) %>% pull(index)]]$fit
dam_ace_fit   <- iam_results[[iam_results_diagnostic %>% filter(model_name == "dam_ace" ) %>%slice(1) %>% pull(index)]]$fit
iam_ae_fit   <- iam_results[[iam_results_diagnostic %>% filter(model_name == "iam_ae" )%>% slice(1) %>% pull(index)]]$fit

mxCompare(qiam_cate_fit, qiam_ace_fit)
mxCompare(qiam_ace_fit, qdam_ace_fit)
mxCompare(qiam_ace_fit, iam_ace_fit)
mxCompare(iam_ace_fit, dam_ace_fit)
mxCompare(iam_ace_fit, iam_ae_fit)

xdam_ade_fit$est_var # biologically plausible

# model of interest (re-run to obtain confidence intervals)
iam_results_diagnostic %>% filter(model_name == "iam_ace")
iam_results_diagnostic %>% filter(model_name == "dam_ace")

iam_ace_fit <- mxRun(iam_ace_mod, intervals = T)
dam_ace_fit  <- mxRun(dam_ace_mod, intervals = T)

iam_ace_est <- summary(iam_ace_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(iam_ace_fit$est_path$result), colnames(iam_ace_fit$est_var$result), colnames(iam_ace_fit$est_cov$result)), model = "iAM-ACE")
dam_ace_est <- summary(dam_ace_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(dam_ace_fit$est_path$result), colnames(dam_ace_fit$est_var$result), colnames(dam_ace_fit$est_cov$result)), model = "dAM-ACE")

est <- rbind(iam_ace_est, dam_ace_est)

# Save table 4
write_csv(est %>% filter(!parameter %in% c("varP_N", "varS_N", "n1s", "n1", "rN_Partner", "rN_FS")) %>% select(model, parameter, estimate, lbound, ubound) %>% mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5)), sprintf("%s/SI/SI.04_table_est_EA.csv", wd_oa))

# SAVE ####
# final models
xam_cdate_mxcompare <- rbind(
  # comparison with am DATE and sibling-shared preferences
  as.data.frame(mxCompare(esat_fit, qiam_cate_fit)),
  as.data.frame(mxCompare(qiam_cate_fit, qiam_ace_fit))[-1,],
  as.data.frame(mxCompare(qiam_ace_fit,  qdam_ace_fit))[-1,],
  as.data.frame(mxCompare(qiam_ace_fit, iam_ace_fit))[-1,],
  as.data.frame(mxCompare(iam_ace_fit, iam_ae_fit))[-1,]
)

xam_cdate_mxcompare <- xam_cdate_mxcompare %>%
  dplyr::select(-c(fitUnits, fit, diffFit, chisq, SBchisq)) %>%
  mutate(minus2LL = round(minus2LL, 0), AIC = round(AIC, 0), diffLL = round(diffLL, 2), p = round(p, 5))

write_csv(xam_cdate_mxcompare, sprintf("%s/SI/SI.04_table_mod_comp_EA.csv", wd_oa))

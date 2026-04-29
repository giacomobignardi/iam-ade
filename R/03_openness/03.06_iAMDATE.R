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
dat_fam2g <- read.csv(sprintf("%s/%s/03.02_dat_estd_openness.csv", wd_noa, wd_noa_data))
dat_fam2g <- dat_fam2g %>% rename(value_tw1 = tw1, value_tw2 = tw2, value_s1 = s1, value_stw1 = stw1, value_stw2 = stw2)

# make data long
dat_fam2g_long <- dat_fam2g %>% pivot_longer(cols = age_tw1:value_stw2, names_to = c(".value", "rel"), names_sep = "_")

# inspect tails (extreme cases may need to be removed)
ifelse(is_cluster == F, hist(dat_fam2g_long$value), print("no GUI in current env"))
nrow(dat_fam2g_long)
# 44725

# TRANSFORMATION & RESIDUALIZATION ####
# residualise for sex and age
dat_fam2g_long$res <- scale(residuals(lm(value ~ sex + age, data = dat_fam2g_long, na.action = na.exclude)))

# pivot wider
dat_fam2g <- dat_fam2g_long %>%
  dplyr::select(fid, twzyg, rel, res) %>%
  pivot_wider(names_from = rel, values_from = res)
dat_fam2g <- dat_fam2g %>% rename(value_tw1 = tw1, value_tw2 = tw2, value_s1 = s1, value_stw1 = stw1, value_stw2 = stw2)

# load saturated model for later model comparison
load(sprintf("%s/%s/03.03_saturated_openness.Rdata", wd_noa, wd_noa_data))
load(sprintf("%s/%s/03.03_saturated_reduced_openness.Rdata", wd_noa, wd_noa_data))

# FUNCTIONS & DATA PREPARATION ####
# source ctd and estd inspired model
source(sprintf("%s/%s/function/iam.date.fun.new.R", wd_oa, wd_oa_scripts))
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
mean_MZ <- c("mean", "mean_TW", "mean", "mean_TW", "mean")
mean_DZ <- c("mean", "mean_TW", "mean", "mean_TW", "mean")

# ADE ####
# first fit classical twin design
# with dominance
ade_mod <- acde.fun(mz, dz, c("value_tw1", "value_tw2"), dp.free = T)
ade_fit <- mxRun(ade_mod, intervals = T)
ade_fit$US$result

# without shared environmental effects
ae_mod <- acde.fun(mz, dz, c("value_tw1", "value_tw2"))
ae_fit <- mxRun(ae_mod, intervals = T)
ae_fit$US$result

# final models
ade_mxcompare <- as.data.frame(mxCompare(ade_fit, ae_fit)) %>%
  dplyr::select(-c(fitUnits, fit, diffFit, chisq, SBchisq)) %>%
  mutate(minus2LL = round(minus2LL, 0), AIC = round(AIC, 0), diffLL = round(diffLL, 2), p = round(p, 5))

ae_est <- summary(ae_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c("V","A","D","C","E"), model = "AE")

# Save table CTD
write_csv(ae_est %>% select(model, parameter, estimate, lbound, ubound) %>% mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5)), sprintf("%s/SI/03.06_tableCTD_est_openness.csv", wd_oa))

# HELPER FUNCTION ####
# Note that sourcing this will overwrite some model specification (e.g. ACE/ADE)
source(sprintf("%s/%s/function/iam.date.fun.helper.R", wd_oa, wd_oa_scripts))

# iAM-DATE ####
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
iam_date_fit <- iam_results[[iam_results_diagnostic %>% filter(model_name == "iam_date" ) %>% pull(index)]]$fit
iam_ade_fit <- iam_results[[iam_results_diagnostic %>% filter(model_name == "iam_ade" ) %>% pull(index)]]$fit
mxCompare(iam_date_fit, iam_ade_fit)

dam_ade_fit <- iam_results[[iam_results_diagnostic %>% filter(model_name == "dam_ade" )%>% slice(1) %>% pull(index)]]$fit
mxCompare(iam_ade_fit, dam_ade_fit)

dam_ae_fit <- iam_results[[iam_results_diagnostic %>% filter(model_name == "dam_ae" ) %>% slice(1) %>% pull(index)]]$fit
mxCompare(dam_ade_fit, dam_ae_fit)

dam_ade_fit$est_var 
dam_ade_fit$est_path
# biologically plausible

# # inspect fit of the model if needed
# mxGetExpected(dam_ade_fit, "covariance")$MZ-cov(mz, use = "pairwise.complete.obs")
# mxGetExpected(dam_ade_fit, "covariance")$DZ-cov(dz, use = "pairwise.complete.obs")
# mxGetExpected(dam_ade_fit, "covariance")$MZ-mxGetExpected(esat_law_fit, "covariance")$mz
# mxGetExpected(dam_ade_fit, "covariance")$DZ-mxGetExpected(esat_law_fit, "covariance")$dz

# model of interest (re-run to obtain confidence intervals)
iam_results_diagnostic %>% filter(model_name == "iam_ade")
iam_results_diagnostic %>% filter(model_name == "dam_ade")

iam_ade_fit <- mxRun(iam_ade_mod, intervals = T)
dam_ade_fit  <- mxRun(dam_ade_mod, intervals = T)

iam_ade_est <- summary(iam_ade_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(iam_ade_fit$est_path$result), colnames(iam_ade_fit$est_var$result), colnames(iam_ade_fit$est_cov$result)), model = "iAM-ADE")
dam_ade_est <- summary(dam_ade_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(dam_ade_fit$est_path$result), colnames(dam_ade_fit$est_var$result), colnames(dam_ade_fit$est_cov$result)), model = "dAM-ADE")

est <- rbind(iam_ade_est, dam_ade_est)

# Save table 4
write_csv(est %>% filter(!parameter %in% c("varP_N", "varS_N", "n1s", "n1", "rN_Partner", "rN_FS")) %>% select(model, parameter, estimate, lbound, ubound) %>% mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5)), sprintf("%s/SI/03.06_table_est_openness.csv", wd_oa))

# SAVE ####
# final models
xam_cdate_mxcompare <- rbind(
  # comparison with am DATE and sibling-shared preferences
  as.data.frame(mxCompare(esat_fit,iam_date_fit)),
  as.data.frame(mxCompare(iam_date_fit, iam_ade_fit))[-1,],
  as.data.frame(mxCompare(iam_ade_fit, dam_ade_fit))[-1,],
  as.data.frame(mxCompare(dam_ade_fit, dam_ae_fit))[-1,]
)

xam_cdate_mxcompare <- xam_cdate_mxcompare %>%
  dplyr::select(-c(fitUnits, fit, diffFit, chisq, SBchisq)) %>%
  mutate(minus2LL = round(minus2LL, 0), AIC = round(AIC, 0), diffLL = round(diffLL, 2), p = round(p, 5))

# save for supplementary table
write_csv(xam_cdate_mxcompare, sprintf("%s/SI/03.06_table_mod_comp_openness.csv", wd_oa))

## rev1.R1.1 including (fixed) measurement error/intra-individual variance ####
# approximate estimated test-retest reliability
load(sprintf("%s/%s/03.01_rYY_openness.Rdata", wd_noa, wd_noa_data))

# Final model including measurement error
dam_ader_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = F,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              d1.free = T, c1.free = F, t1.free = F,
                              d1s.free = T, c1s.free = F, t1s.free = F,
                              d1.val = .6, a1.val = .6, c1.val = .0, t1.val = .0, e1.val = rYY$ICC,
                              # test/retest
                              r1.free = TRUE, 
                              ryy = rYY$ICC,
                              name = "dam_ader"
)

dam_ader_fit <- mxRun(dam_ader_mod, intervals = T)
dam_ader_fit$est_var
dam_ader_est <- summary(dam_ader_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(dam_ader_fit$est_path$result), colnames(dam_ader_fit$est_var$result), colnames(dam_ader_fit$est_cov$result)), model = "dAM-ADE")

# Save table 5
write_csv(dam_ader_est %>% filter(!parameter %in% c("varP_N", "varS_N", "n1s", "n1", "rN_Partner", "rN_FS")) %>% select(model, parameter, estimate, lbound, ubound) %>% mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5)), sprintf("%s/SI/03.06_table_est_error_openness.csv", wd_oa))

## rev1.RX.X 0th generation of assortment ####
# Final model assuming first generation of assortment
am_ade_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = F,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = T, c1.free = F, t1.free = F,
                            d1s.free = T, c1s.free = F, t1s.free = F,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            # generation of assortment
                            U.f = 0,
                            name = "am_ade"
)

am_ade_fit <- mxRun(am_ade_mod, intervals = T)
am_ade_fit$est_var
am_ade_est <- summary(am_ade_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(am_ade_fit$est_path$result), colnames(am_ade_fit$est_var$result), colnames(am_ade_fit$est_cov$result)), model = "ADE (0th gen am)")

# Save table 6
write_csv(am_ade_est %>% filter(!parameter %in% c("varP_N", "varS_N", "n1s", "n1", "rN_Partner", "rN_FS")) %>% select(model, parameter, estimate, lbound, ubound) %>% mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5)), sprintf("%s/SI/03.06_table_est_0tham_openness.csv", wd_oa))

## rev1.R2.8 dominance under random mating ####
iam_results_diagnostic %>% filter(model_name == "ade")
ade_fit  <- mxRun(ade_mod, intervals = T)
ade_est <- summary(ade_fit)$CI %>%
  as.data.frame() %>%
  mutate(parameter = c(colnames(ade_fit$est_path$result), colnames(ade_fit$est_var$result), colnames(ade_fit$est_cov$result)), model = "ADE")

# Save table 7
write_csv(ade_est %>% filter(!parameter %in% c("varP_N", "varS_N", "n1s", "n1", "rN_Partner", "rN_FS")) %>% select(model, parameter, estimate, lbound, ubound) %>% mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5)), sprintf("%s/SI/03.06_table_est_noam_openness.csv", wd_oa))

# EPISTASIS####
# iam_nate with 2 to 10 AxA 
# fit the model with direct assortment and dominance but not twin-specific effects
dam_ane_est <- c()
for (i in seq(2, 10, 1)) {
  dam_ane_mod <- iam.date.fun(mz,
                               dz,
                               sel_vars,
                               indirect.assortment = F,
                               # unconstrained means
                               mean.MZ = mean_MZ,
                               mean.DZ = mean_DZ,
                               # focal paths
                               d1.free = F,
                               c1.free = F,
                               t1.free = F,
                               n1.free = T,
                               # sorting path
                               d1s.free = F,
                               c1s.free = F,
                               t1s.free = F,
                               n1s.free = T,
                               # starting values
                               d1.val = .0,
                               a1.val = .8,
                               c1.val = .0,
                               t1.val = .0,
                               e1.val = .3,
                               # AxA
                               ni = i,
                               # name of the model
                               name = "xdam_ane"
  )
  dam_ane_fit <- mxTryHard(dam_ane_mod, intervals = T)
  dam_ane_est_i <- summary(dam_ane_fit)$CI %>%
    as.data.frame() %>%
    mutate(parameter = c(colnames(dam_ane_fit$est_path$result), colnames(dam_ane_fit$est_var$result), colnames(dam_ane_fit$est_cov$result)), model = "dAM-AIE", ri = i)
  dam_ane_est <- rbind(dam_ane_est, dam_ane_est_i)
}

# save for later
nate_est_tb <- dam_ane_est %>%
  filter(parameter %in% c("varP_A", "varP_N", "varP_E", "rN_Partner", "rN_FS")) %>%
  mutate(
    parameter = case_match(
      parameter,
      "varP_E" ~ "E",
      "varP_N" ~ "N",
      "varP_A" ~ "A",
      "rN_Partner" ~ "rN_Partner",
      "rN_FS" ~ "rN_FS"
    ),
    parameter = factor(parameter, levels = c("A", "N", "E", "rN_Partner","rN_FS"))
  ) %>%
  rename(axa = ri) %>%
  mutate(rNA = 1/(2^axa)) %>%
  select(model, axa, rNA, parameter, estimate, lbound, ubound) %>%
  arrange(axa) %>%
  mutate(lbound = round(lbound, 5), estimate = round(estimate, 5), ubound = round(ubound, 5))

# Save table 8
write_csv(nate_est_tb, sprintf("%s/SI/03.06_table_est_nonadditive_openness.csv", wd_oa))
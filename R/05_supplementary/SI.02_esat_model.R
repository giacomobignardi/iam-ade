# SCRIPT:        esat_model.R
# AUTHOR:        G.B. (Giaco)
# ADAPTED FROM:  https://hermine-maes.squarespace.com/#/one/ saturated models
# PURPOSE:       Fit saturated and semi-constrained saturated models

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
# load extended twin family data
dat_fam2g <- read.csv(sprintf("%s/%s/SI.01_dat_estd_EA.csv", wd_noa, wd_noa_data))
dat_fam2g <- dat_fam2g %>% rename(value_tw1 = tw1, value_tw2 = tw2, value_s1 = s1, value_stw1 = stw1, value_stw2 = stw2)

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
  )) %>%   # scale within sex (for consistency Sunde et al. 2025 https://doi.org/10.1038/s41467-025-60483-0)
  group_by(sex) %>% 
  mutate(z = as.numeric(scale(value))) %>% 
  ungroup() 

# dat_fam2g_long$res <- scale(residuals(lm(value ~ sex + age, data = dat_fam2g_long, na.action = na.exclude)))
# ifelse(is_cluster == F, hist(dat_fam2g_long$res), print("no GUI in current env"))
# dat_fam2g_long$resrank <- qnorm((rank(dat_fam2g_long$res,na.last="keep")-0.5)/sum(!is.na(dat_fam2g_long$res)))
# dat_fam2g_long$resrankres <- scale(residuals(lm(resrank ~ sex + age, data = dat_fam2g_long, na.action = na.exclude)))[,1]
# ifelse(is_cluster == F, hist(dat_fam2g_long$resrankres), print("no GUI in current env"))

# pivot wider
dat_fam2g <- dat_fam2g_long %>%
  dplyr::select(fid, twzyg, rel, z) %>%
  pivot_wider(names_from = rel, values_from = z)
dat_fam2g <- dat_fam2g %>% rename(value_tw1 = tw1, value_tw2 = tw2, value_s1 = s1, value_stw1 = stw1, value_stw2 = stw2)

# select Variables for Analysis
sel_vars <- paste0("value", "_", c("stw1", "tw1", "s1", "tw2", "stw2"))
# select Data for Analysis
mz <- subset(dat_fam2g, twzyg == "MZ", sel_vars)
dz <- subset(dat_fam2g, twzyg == "DZ", sel_vars)

# observed correlations
cov(mz, use = "pairwise.complete.obs")
cov(dz, use = "pairwise.complete.obs")

psych::pairwiseCount(mz)
psych::pairwiseCount(dz)
# DYADS ####
# number of complete observation
mz_sample <- psych::pairwiseCount(mz)
dz_sample <- psych::pairwiseCount(dz)

colnames(mz_sample) <- sel_vars
rownames(mz_sample) <- sel_vars
colnames(dz_sample) <- sel_vars
rownames(dz_sample) <- sel_vars

# cor MZ twin
nmz_comp_pair <- mz_sample["value_tw1", "value_tw2"]

# cor DZ twin
ndz_comp_pair <- dz_sample["value_tw1", "value_tw2"]

# cor FS
nfs_comp_pair <- mz_sample["value_s1", "value_tw1"] + mz_sample["value_s1", "value_tw2"] + dz_sample["value_s1", "value_tw1"] + dz_sample["value_s1", "value_tw2"]

# cor spouses (am)
nstw_comp_pair <- mz_sample["value_tw1", "value_stw1"] + dz_sample["value_tw2", "value_stw2"] + mz_sample["value_tw2", "value_stw2"] + dz_sample["value_tw1", "value_stw1"]

# cor MZ twin in law
nmzil_comp_pair <- mz_sample["value_tw1", "value_stw2"] + mz_sample["value_tw2", "value_stw1"]

# cor DZ twin in law
ndzil_comp_pair <- dz_sample["value_tw1", "value_stw2"] + dz_sample["value_tw2", "value_stw1"]
nfsil_comp_pair <- mz_sample["value_s1", "value_stw1"] + mz_sample["value_s1", "value_stw2"] + dz_sample["value_s1", "value_stw1"] + dz_sample["value_s1", "value_stw2"]

# cor twin in co-law
nmzcil_comp_pair <- mz_sample["value_stw1", "value_stw2"]

# cor twin in co-law
ndzcil_comp_pair <- dz_sample["value_stw1", "value_stw2"]

# save for later
n_compair <- data.frame(
  # cor MZ 
  nmz_comp_pair,
  # cor DZ 
  ndz_comp_pair,
  # cor FS 
  nfs_comp_pair,
  # cor spouses (am)
  nstw_comp_pair,
  # cor MZ twin in law
  nmzil_comp_pair,
  # cor DZ twin in law
  ndzil_comp_pair,
  # cor FS twin in law
  nfsil_comp_pair,
  # cor mz in co-law
  nmzcil_comp_pair,
  # cor dz in co-law
  ndzcil_comp_pair
)

# save complete pair to add to later table
write_csv(n_compair, sprintf("%s/%s/SI.02_pair_number_EA.csv", wd_noa, wd_noa_data))

# check that we are capturing all relationships (total number of pairwise complete observation)
sum(mz_sample[lower.tri(mz_sample)]) + sum(dz_sample[lower.tri(dz_sample)]) ==

  # cor MZ 
  nmz_comp_pair +
    # cor DZ 
    ndz_comp_pair +
    # cor FS 
    nfs_comp_pair +
    # cor spouses (am)
    nstw_comp_pair +
    # cor MZ twin in law
    nmzil_comp_pair +
    # cor DZ twin in law
    ndzil_comp_pair +
    # cor FS twin in law
    nfsil_comp_pair +
    # cor twin in co-law
    nmzcil_comp_pair +
    ndzcil_comp_pair

# total number of individuals
sum(diag(mz_sample)) + sum(diag(dz_sample))
# 11206

# MODEL ####
# set starting Values
m.va <- 0 # start value for means
v.va <- 1 # start value for variance
lb.Va <- .0001 # lower bound for variance
## Model specification####
# create algebra for expected mean matrices
meanmz <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = m.va, labels = c("mmzstw1", "mmz1", "mmzfs", "mmz2", "mmzstw2"), name = "meanmz")
meandz <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = m.va, labels = c("mdzstw1", "mdz1", "mdzfs", "mdz2", "mdzstw2"), name = "meandz")

# create algebra for expected variance/covariance matrices
covmz <- mxMatrix( type="Full", nrow=5, ncol=5, free=TRUE, values=diag(5),
                   labels=c("vmzstw1",   "cmzstw1", "cmzfsil1",  "cmzil2",   "cmzstwcil",
                            "cmzstw1",   "vmz1",    "cmzfs1",    "cmz",      "cmzil1",
                            "cmzfsil1",  "cmzfs1",  "vmzfs",     "cmzfs2",  "cmzfsil2",
                            "cmzil2",    "cmz",     "cmzfs2",   "vmz2",     "cmzstw2",
                            "cmzstwcil", "cmzil1",  "cmzfsil2",  "cmzstw2",  "vmzstw2"), name="covmz" )

# create data objects for multiple groups
covdz <- mxMatrix( type="Full", nrow=5, ncol=5, free=TRUE, values=diag(5),
                   labels=c("vdzstw1",   "cdzstw1", "cdzfsil1",   "cdzil2",  "cdzstwcil",
                            "cdzstw1",   "vdz1",    "cdzfs1",     "cdz",     "cdzil1",
                            "cdzfsil1",  "cdzfs1",  "vdzfs",      "cdzfs2", "cdzfsil2",
                            "cdzil2",    "cdz",     "cdzfs2",     "vdz2",    "cdzstw2",
                            "cdzstwcil", "cdzil1",  "cdzfsil2",   "cdzstw2", "vdzstw2"), name="covdz" )

# correlations
cormz <- mxAlgebra(cov2cor(covmz), name = "cormz")
cordz <- mxAlgebra(cov2cor(covdz), name = "cordz")

# create data objects for multiple groups
datamz <- mxData(observed = mz, type = "raw")
datadz <- mxData(observed = dz, type = "raw")

# create expectation objects for multiple groups
expmz <- mxExpectationNormal(covariance = "covmz", means = "meanmz", dimnames = sel_vars)
expdz <- mxExpectationNormal(covariance = "covdz", means = "meandz", dimnames = sel_vars)
funML <- mxFitFunctionML()

# create model objects for multiple groups
modelmz <- mxModel(meanmz, covmz, cormz, datamz, expmz, funML, name = "mz")
modeldz <- mxModel(meandz, covdz, cordz, datadz, expdz, funML, name = "dz")

multi <- mxFitFunctionMultigroup(c("mz", "dz"))

# create confidence interval objects
ciCov <- mxCI(c("mz.covmz", "dz.covdz"))
ciCor <- mxCI(c("mz.cormz", "dz.cordz"))
ciMean <- mxCI(c("mz.meanmz", "dz.meandz"))

# build the extended saturated model with confidence intervals
esat_mod <- mxModel("saturated", modelmz, modeldz, multi, ciCov, ciCor, ciMean)

# FIT ####
# run saturated model
esat_fit <- mxRun(esat_mod, intervals = F)
esat_sumy <- summary(esat_fit)
esat_sumy$parameters
esat_sumy$CI
esat_sumy$parameters %>% arrange(name)

save(esat_fit, file = sprintf("%s/%s/SI.02_saturated_EA.Rdata", wd_noa, wd_noa_data))

# CONSTRAINS ####
# constrain mean and variances and covariances to have a reduced number of correlations
# constrain means across twin order
esat_mbo_mod <- mxModel(esat_fit, name = "means twin birth order")

# equate means across twin order
# twin order means (it is know that first-born are somewhat taller then second-born twins)
esat_mbo_mod <- omxSetParameters(esat_mbo_mod, label = c("mmz1", "mmz2"), free = TRUE, values = 1, newlabels = "mmz")
esat_mbo_mod <- omxSetParameters(esat_mbo_mod, label = c("mdz1", "mdz2"), free = TRUE, values = 1, newlabels = "mdz")

# test if we can equate means
esat_mbo_fit <- mxRun(esat_mbo_mod, intervals = F)
mxCompare(esat_fit, esat_mbo_fit) 

# constrain means and variances across  pairs
esat_mvar_mod <- mxModel(esat_mbo_mod, name = "means and variances order")

# equate means and variances to obtain summary statistics on the following:
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("mdzstw1", "mdzstw2"), free = TRUE, values = 1, newlabels = "mdzstw")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("mmzstw1", "mmzstw2"), free = TRUE, values = 1, newlabels = "mmzstw")

# twin order var
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vmz1", "vmz2"), free = TRUE, values = 1, newlabels = "vmz")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vdz1", "vdz2"), free = TRUE, values = 1, newlabels = "vdz")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vdzstw1", "vdzstw2"), free = TRUE, values = 1, newlabels = "vdzstw")
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vmzstw1", "vmzstw2"), free = TRUE, values = 1, newlabels = "vmzstw")

# sibling  means
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("mmzfs", "mdzfs"), free = TRUE, values = 1, newlabels = "mfs")
# sibling  variances
esat_mvar_mod <- omxSetParameters(esat_mvar_mod, label = c("vmzfs", "vdzfs"), free = TRUE, values = 1, newlabels = "vfs")

# fit and test model fit
esat_mvar_fit <- mxRun(esat_mvar_mod, intervals = F)
mxCompare(esat_fit, esat_mvar_fit)

# extract parameters
esat_mvar_sumy <- summary(esat_mvar_fit)
esat_mvar_sumy$parameters

# element to later augment table 1.0
var_mz <- esat_mvar_sumy$parameters %>%
  filter(name == "vmz") %>%
  pull(Estimate)
var_dz <- esat_mvar_sumy$parameters %>%
  filter(name == "vdz") %>%
  pull(Estimate)
var_sib <- esat_mvar_sumy$parameters %>%
  filter(name == "vfs") %>%
  pull(Estimate)
mean_mz <- esat_mvar_sumy$parameters %>%
  filter(name == "mmz") %>%
  pull(Estimate)
mean_dz <- esat_mvar_sumy$parameters %>%
  filter(name == "mdz") %>%
  pull(Estimate)
mean_sib <- esat_mvar_sumy$parameters %>%
  filter(name == "mfs") %>%
  pull(Estimate)
var_mzpartners <- esat_mvar_sumy$parameters %>%
  filter(name == "vmzstw") %>%
  pull(Estimate)
mean_mzpartners <- esat_mvar_sumy$parameters %>%
  filter(name == "mmzstw") %>%
  pull(Estimate)
var_dzpartners <- esat_mvar_sumy$parameters %>%
  filter(name == "vdzstw") %>%
  pull(Estimate)
mean_dzpartners <- esat_mvar_sumy$parameters %>%
  filter(name == "mdzstwZ") %>%
  pull(Estimate)

# equate means and variances to test for various differences
# equate mean across zygosities
esat_mtw1_mod <- mxModel(esat_mvar_mod, name = "mean across zyg")
esat_mtw1_mod <- omxSetParameters(esat_mtw1_mod, label = c("mmz", "mdz"), free = TRUE, values = 1, newlabels = "mtw")
esat_mtw1_fit <- mxRun(esat_mtw1_mod, intervals = F)
mxCompare(esat_fit, esat_mtw1_fit) 


# equate mean across zygosities and siblings
esat_msib1_mod <- mxModel(esat_mtw1_mod, name = "mean across zyg and sib")
esat_msib1_mod <- omxSetParameters(esat_mtw1_mod, label = c("mtw", "mfs"), free = TRUE, values = 1, newlabels = "msib")
esat_msib1_fit <- mxRun(esat_msib1_mod, intervals = F)
mxCompare(esat_fit, esat_msib1_fit) # we cannot constrain means


# equate variances across zygosities
esat_vartw1_mod <- mxModel(esat_mvar_mod, name = "var across zyg and sib")
esat_vartw1_mod <- omxSetParameters(esat_vartw1_mod, label = c("vmz", "vdz", "vfs"), free = TRUE, values = 1, newlabels = "vsib")
esat_vartw1_fit <- mxRun(esat_vartw1_mod, intervals = F)
esat_vartw1_fit <- mxTryHard(esat_vartw1_mod, intervals = F)
mxCompare(esat_fit, esat_vartw1_fit)

# equate mean and variances across zygosities for spouses of twins
esat_mvastw_mod <- mxModel(esat_vartw1_mod, name = "mean and var across zyg for partners")
esat_mvastw_mod <- omxSetParameters(esat_mvastw_mod, label = c("mmzstw", "mdzstw"), free = TRUE, values = 1, newlabels = "mstw")
esat_mvastw_mod <- omxSetParameters(esat_mvastw_mod, label = c("vmzstw", "vdzstw"), free = TRUE, values = 1, newlabels = "vstw")
esat_mvastw_fit <- mxRun(esat_mvastw_mod, intervals = F)
mxCompare(esat_fit, esat_mvastw_fit)

# equate variances
esat_var1_mod <- mxModel(esat_mvastw_mod, name = "var partners and sib")
esat_var1_mod <- omxSetParameters(esat_var1_mod, label = c("vstw", "vsib"), free = TRUE, values = 1, newlabels = "vp")
esat_var1_fit <- mxRun(esat_var1_mod, intervals = F)
esat_var1_fit <- mxTryHard(esat_var1_fit, intervals = F)
mxCompare(esat_fit, esat_var1_fit)  # we cannot constrain variances

# extract parameters
esat_mvastw_sumy <- summary(esat_mvastw_fit)
esat_mvastw_sumy$parameters

# equate means
esat_mean1_mod <- mxModel(esat_mvastw_mod, name = "mean partners and sib")
esat_mean1_mod <- omxSetParameters(esat_mean1_mod, label = c("mstw", "mfs"), free = TRUE, values = 1, newlabels = "mp")
esat_mean1_fit <- mxRun(esat_mean1_mod, intervals = F)
esat_mean1_fit <- mxTryHard(esat_mean1_fit, intervals = F)
mxCompare(esat_fit, esat_mean1_fit) # we cannot constrain means

summary(esat_mean1_fit)

# sibling correlations
esat_sib_mod <- mxModel(esat_mvastw_mod, name = "covar full siblings across zyg")
esat_sib_mod <- omxSetParameters(esat_sib_mod, label = c("cdzfs1", "cdzfs2", "cmzfs1", "cmzfs2"), free = TRUE, values = .1, newlabels = "cfs")
esat_sib_fit <- mxRun(esat_sib_mod, intervals = F)
esat_sib_fit <- mxTryHard(esat_sib_mod, intervals = F)
mxCompare(esat_fit, esat_sib_fit)

# assortative mating across zygosities
esat_asmtw_mod <- mxModel(esat_sib_mod, name = "cov partners across zyg")
esat_asmtw_mod <- omxSetParameters(esat_asmtw_mod, label = c("cmzstw1", "cmzstw2", "cdzstw1", "cdzstw2"), free = TRUE, values = .1, newlabels = "cprtw")
esat_asmtw_fit <- mxRun(esat_asmtw_mod, intervals = F)
esat_asmtw_fit <- mxTryHard(esat_asmtw_mod, intervals = F)
mxCompare(esat_fit, esat_asmtw_fit)

# in law
esat_law_mod <- mxModel(esat_asmtw_mod, name = "cov in law across zyg")
esat_law_mod <- omxSetParameters(esat_law_mod, label = c("cmzil1", "cmzil2"), free = TRUE, values = .1, newlabels = "cmzinlaw")
esat_law_mod <- omxSetParameters(esat_law_mod, label = c("cdzil1", "cdzil2"), free = TRUE, values = .1, newlabels = "cdzinlaw")
esat_law_mod <- omxSetParameters(esat_law_mod, label = c("cdzfsil1", "cdzfsil2", "cmzfsil1", "cmzfsil2"), free = TRUE, values = .1, newlabels = "cfsinlaw")
esat_law_fit <- mxRun(esat_law_mod, intervals = F)
esat_law_fit <- mxTryHard(esat_law_mod, intervals = F)
mxCompare(esat_fit, esat_law_fit)
esat_law_sumy <- summary(esat_law_fit)
esat_law_sumy$parameters %>% arrange(name)


# TEST ####
# Sibs
esat_csib1_mod <- mxModel(esat_law_mod, name = "test1: full-sibling vs dz")
esat_csib1_mod <- omxSetParameters(esat_csib1_mod, label = c("cdz", "cfs"), free = TRUE, values = .1, newlabels = "cfs")
esat_csib1_fit <- mxRun(esat_csib1_mod, intervals = F)
# esat_csib1_fit <- mxTryHard(esat_csib1_mod, intervals = F)
mxCompare(esat_law_fit, esat_csib1_fit)

# AM (DZ-FS in law)
esat_am2_mod <- mxModel(esat_law_mod, name = "test2: dz vs fs in law")
esat_am2_mod <- omxSetParameters(esat_am2_mod, label = c("cdzinlaw", "cfsinlaw"), free = TRUE, values = .1, newlabels = "csibinlaw")
esat_am2_fit <- mxRun(esat_am2_mod, intervals = T)
esat_am2_fit <- mxTryHard(esat_am2_fit, intervals = T)
mxCompare(esat_law_fit, esat_am2_fit)

# SAVE ####
save(esat_am2_fit, file = sprintf("%s/%s/SI.02_saturated_reduced_EA.Rdata", wd_noa, wd_noa_data))

## mx compare list ####
esat_model_comparsion <- mxCompare(esat_fit, c(
  esat_mbo_fit,
  esat_mvar_fit,
  esat_mtw1_fit,
  esat_vartw1_fit,
  esat_mvastw_fit,
  esat_var1_fit,
  esat_mean1_fit,
  esat_sib_fit,
  esat_asmtw_fit,
  esat_law_fit,
  # test sib differences
  esat_csib1_fit,
  # test sib in-law differences
  esat_am2_fit
))

esat_model_comparsion <- esat_model_comparsion %>%
  as.data.frame() %>%
  dplyr::select(-c(fitUnits, fit, diffFit, chisq, SBchisq)) %>%
  mutate(minus2LL = round(minus2LL, 0), AIC = round(AIC, 0), diffLL = round(diffLL, 2), p = round(p, 5))

## supplementary table ####
write_csv(as.data.frame(esat_model_comparsion), sprintf("%s/SI/SI.03_table_esat_mxcompare_EA.csv", wd_oa))

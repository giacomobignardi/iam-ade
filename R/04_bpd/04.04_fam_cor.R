# SCRIPT:  fam_cor.R
# AUTHOR:  G.B. (Giaco)
# PURPOSE: Extract family sample size and correlations between family members

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
# load estimates from the constrained sat models
load(file = sprintf("%s/%s/04.03_saturated_reduced_bpd.Rdata", wd_noa, wd_noa_data))

# load number of complete pairs
n_compair <- read.csv(sprintf("%s/%s/04.03_pair_number_bpd.csv", wd_noa, wd_noa_data))

# FAMILY COV ####
esat_am2_sumy <- summary(esat_am2_fit)

fam_cov <- esat_am2_sumy$parameters %>%
  filter(grepl("^c", name)) %>%
  dplyr::select(name, Estimate, Std.Error) %>%
  mutate(name = recode(name,
                       "cprtw"      = "am",
                       "cmz"        = "mz",
                       "cdz"        = "dz",
                       "cfs"        = "sib",
                       "cmzinlaw"   = "mz-in-law",
                       "csibinlaw"  = "sib-in-law",
                       "cmzstwcil"  = "mz-co-in-law",
                       "cdzstwcil"  = "dz-co-in-law"
  ))

# FAMILY COR (CIs)####
fam_cor <- esat_am2_sumy$CI %>%
  rownames_to_column() %>%
  mutate(type = recode(rowname,
                       "mz.cormz[2,1]" = "ram",           # Partners
                       "mz.cormz[4,2]" = "rmz",           # MZ twins
                       "dz.cordz[4,2]" = "rdz",           # DZ twins
                       "mz.cormz[3,2]" = "rsib",          # Full siblings
                       "mz.cormz[4,1]" = "rmz-in-law",    # MZ in-law
                       "dz.cordz[3,1]" = "rsib-in-law",   # Siblings in-law
                       "mz.cormz[5,1]" = "rmz-co-in-law", # MZ co-in-law
                       "dz.cordz[5,1]" = "rdz-co-in-law"  # DZ co-in-law
  )) %>%
  filter(type %in% c(
    "ram", "rmz", "rdz", "rsib",
    "rmz-in-law", "rsib-in-law", "rmz-co-in-law", "rdz-co-in-law"
  ))

# SAMPLE SIZE ####
# Combine DZ and full-sibling in-law pairs
n_compair <- n_compair %>%
  mutate(nsil_comp_pair = ndzil_comp_pair + nfsil_comp_pair)

n_compair_long <- n_compair %>%
  rename(
    "ram"            = nstw_comp_pair,
    "rmz"            = nmz_comp_pair,
    "rdz"            = ndz_comp_pair,
    "rsib"           = nfs_comp_pair,
    "rmz-in-law"     = nmzil_comp_pair,
    "rsib-in-law"    = nsil_comp_pair,
    "rmz-co-in-law"  = nmzcil_comp_pair,
    "rdz-co-in-law"  = ndzcil_comp_pair
  ) %>%
  dplyr::select(-c(ndzil_comp_pair, nfsil_comp_pair)) %>%
  pivot_longer(
    cols      = rmz:"rsib-in-law",
    names_to  = "type",
    values_to = "n"
  )

## Merge correlations and sample sizes ####
fam_cor_fin <- merge(fam_cor, n_compair_long, by = "type", all = TRUE) %>%
  mutate(
    name = substr(type, 2, nchar(type)),
    type = factor(type, levels = rev(c(
      "rmz", "rdz", "rsib",
      "rmz-in-law", "rsib-in-law", "rmz-co-in-law", "rdz-co-in-law",
      "ram"
    )))
  ) %>% dplyr::select(-c(rowname, note))

## Combine covariances and correlations into final table ####
fam_est <- merge(fam_cov, fam_cor_fin, by = "name") %>%
  as.data.frame() %>%
  mutate(name = sapply(strsplit(name, "\\."), `[`, 1)) %>%
  arrange(desc(type)) %>%
  dplyr::select(name, n, Estimate, Std.Error, estimate, lbound, ubound) %>%
  rename(
    relationship = name,
    cov          = Estimate,
    SE           = Std.Error,
    r            = estimate
  ) %>%
  mutate(
    across(c(cov, SE, r, lbound, ubound), ~ round(.x, 3)),
    relationship = recode(relationship,
                          "am"           = "Spouses",
                          "mz"           = "MZ",
                          "dz"           = "DZ",
                          "sib"          = "Full sibling",
                          "mz-in-law"    = "MZ in-law",
                          "sib-in-law"   = "Full sibling in-law",
                          "mz-co-in-law" = "MZ co-in-law",
                          "dz-co-in-law" = "DZ co-in-law"
    )
  )

# save
write_csv(fam_cor_fin, sprintf("%s/%s/04.04_fam_cor_bpd.csv", wd_noa, wd_noa_data))
write_csv(fam_est, sprintf("%s/SI/04.04_table_covr_bpd.csv", wd_oa))

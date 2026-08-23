# SCRIPT:          SI_table.R
# AUTHOR:          G.B. (Giaco)
# ADAPTED FROM:    Bignardi et al. (2025, Nat Commun) https://github.com/giacomobignardi/h2_BMRQ/blob/main/R/04_r1_MRS_facets_correlations.R 
# PURPOSE:         Prepare the Supplementary File

rm(list = ls())

library(tidyverse) # Data wrangling
library(openxlsx)  # Prepare excel file

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

# OUTPUT ####
sit <- data.frame(number = c("02", "03", "04"), trait  = c("height", "openness", "bpd"))

st_list <- c()
for ( i in 1:nrow(sit)){
st1 <- read.csv(sprintf("%s/SI/%s.02_table_dem_%s.csv", wd_oa, sit[i,1], sit[i,2]))
st2 <- read.csv(sprintf("%s/SI/%s.04_table_covr_%s.csv", wd_oa, sit[i,1], sit[i,2]))
st3 <- read.csv(sprintf("%s/SI/08_table_implied_cor_%s.csv", wd_oa, sit[i,2]))
st4 <- read.csv(sprintf("%s/SI/%s.03_table_esat_mxcompare_%s.csv", wd_oa, sit[i,1], sit[i,2]))
st5 <- read.csv(sprintf("%s/SI/%s.05_table_lm_conv_%s.csv", wd_oa, sit[i,1], sit[i,2]))
st6 <- read.csv(sprintf("%s/SI/%s.06_table_mod_comp_%s.csv", wd_oa, sit[i,1], sit[i,2]))
st7 <- read.csv(sprintf("%s/SI/%s.06_table_est_%s.csv", wd_oa, sit[i,1], sit[i,2]))
st8 <- read.csv(sprintf("%s/SI/%s.06_table_est_0tham_%s.csv", wd_oa, sit[i,1], sit[i,2]))
st9 <- read.csv(sprintf("%s/SI/%s.06_table_est_error_%s.csv", wd_oa, sit[i,1], sit[i,2]))
st10 <- read.csv(sprintf("%s/SI/%s.06_table_est_nonadditive_%s.csv", wd_oa, sit[i,1], sit[i,2]))
# order of SI st1, st5, st2, st3, st4, st6
st_list[[sit[i,2]]]  <- list(st1, st2, st3, st4, st5, st6, st7, st8, st9, st10)
}

st_list_df <- list()
for (i in 1:10) {
  df <- NULL
  for (j in 1:length(st_list)) {
    if (is.null(df)) {
      df <- st_list[[j]][[i]] %>% mutate(trait = names(st_list)[j])
    } else {
      df <- bind_rows(df, st_list[[j]][[i]] %>% mutate(trait = names(st_list)[j]))
    }
  }
  st_list_df[[i]] <- df
}

# create a reference table
st0 <- data.frame(sheet = c(paste0("S", seq(0, 10, 1))),
                  content = c(
                    "Overview of sheet contents: provides information about the contents of all other sheets",
                    "Complete data per family member: displays the number of individuals with complete data for each family member",
                    "Observed correlations between family members: displays correlations (and covariances) across family relationship types",
                    "Implied correlations between family members: displays expected correlations across family relationship types as implied by a model with direct assortative mating",
                    "Model comparison — saturated model: contains results from the comparison against the saturated model",
                    "Estimates for convergence: contains results from the linear model with relationship-length × partner interactions",
                    "Model comparison — iAM-ADE model: contains results from the comparison against the iAM-ADE model",
                    "Parameter estimates — iAM-ADE model: contains parameter estimates from the most parsimonious iAM-ADE model",
                    "Parameter estimates — iAM-ADE model (0th generation of assortment): contains parameter estimates from the most parsimonious iAM-ADE model for the 0th generation of assortment",
                    "Parameter estimates — iAM-ADE model with fixed measurement error: contains parameter estimates from the most parsimonious iAM-ADE model disattenuated for measurement error",
                    "Alternative estimates — dAM-ANE across rNA values: presents alternative parameter estimates from the dAM-ANE model across varying rNA values"
                  ))
st_list_df <- c(list(st0), st_list_df)

# create a workbook
wb <- createWorkbook()

# create a list for sheet names
sn <- c("0.overview", 
        "1.n", 
        "2.obsv_cor", 
        "3.impl_cor", 
        "4.sat_mod_comp", 
        "5.conv_test", 
        "6.mod_comp", 
        "7.mod_estm", 
        "8.mod_estm_0tham", 
        "9.mod_estm_meas_err",
        "10.mod_estm_epistasis")

# add each correlation matrix as a sheet
for (i in 1:length(st_list_df)) {
  sheet_name <- sn[i]
  addWorksheet(wb, sheet_name)  # Add a new sheet
  writeData(wb, sheet = sheet_name, x = st_list_df[[i]]) 
}

# add note to sup 5
note_sup4 <- "Notes. Bold entries indicate the final semi-constrained model. Where mean constraints resulted in model deterioration, subsequent models freely estimated means. Despite deterioration in fit for some intermediate constrainted models, comparisons between the fully saturated and final semi-constrained models indicate no overall deterioration in fit of the final model."
note_sup5 <- "Notes. Interaction estimates are derived from the linear model on scaled partner trait values (Par1, Par2) and self-reported relationship duration (Duration) : Par1 = α + β₁Par2 + β₂Duration + β₃(Par2 × Duration) + ε; estimates are for β₃."
note_sup6 <- "Notes. For BPD features, the xiam_ade model returned an OpenMx status code 6 (red) and was therefore excluded from model comparisons; the xiam_date model returned status code 1 (green) and was retained. For further information on OpenMx status codes, see: https://openmx.ssri.psu.edu/wiki/errors"
note_sup7 <- "Notes. The variances (var) are standardised, with the exception of varP = variance of the focal phenotype and varS = variance of the sorting factor; the parameter labels for the path coefficients to the sorting factor end in 's' (e.g., as); models starting with 'x' include cross-trait assortative mating on sibling-shared preferences; iAM = indirect assortative mating; dAM = direct assortative mating."
note_sup8 <- "Notes. As in Supplementary Table 7, but with parameter estimates assuming the population is in the first generation of assortment."
note_sup9 <- "Notes. Parameters including the term '_bw_' are obtained by removing within-individual variance from the overall variance. The xdam_ade model for BPD features returned status code 1 (green) and was retained; MZ = Monozygotic; DZ = Dizygotic; FS = Full Siblings; In = In-laws; CIn = Co-In-laws; cor_Partner_P = correlation between spouses on the focal phenotype; cor_Partner_S = correlation between spouses on the sorting factor (i.e., corrected for within-individual differences)."
note_sup10 <- "Notes. Parameters under different orders of epistasis (from AxA to AxAxAxAxAxAxAxAxAxA). 
rNA = uncorrected interaction-deviation correlation between full siblings; 
rNA_FS = implied interaction-deviation correlation between full siblings under assortment; 
rNA_Partner = implied interaction-deviation correlation between partners;
N = non-additive (interaction-deviation) genetic variance component."

writeData(wb, "4.sat_mod_comp",  note_sup4, startCol = 1, startRow = 5 + nrow(st_list_df[[4+1]]))
writeData(wb, "5.conv_test",  note_sup5, startCol = 1, startRow = 5 + nrow(st_list_df[[5+1]]))
writeData(wb, "6.mod_comp",  note_sup6, startCol = 1, startRow = 5 + nrow(st_list_df[[6+1]]))
writeData(wb, "7.mod_estm", note_sup7, startCol = 1, startRow = 5 + nrow(st_list_df[[7+1]]))
writeData(wb, "8.mod_estm_0tham", note_sup8, startCol = 1, startRow = 5 + nrow(st_list_df[[8+1]]))
writeData(wb, "9.mod_estm_meas_err", note_sup9, startCol = 1, startRow = 5 + nrow(st_list_df[[9+1]]))
writeData(wb, "10.mod_estm_epistasis", note_sup10, startCol = 1, startRow = 5 + nrow(st_list_df[[10+1]]))

# save Supplementary file
saveWorkbook(wb,sprintf("%s/SI/Supplementary_Table.xlsx", wd_oa), overwrite = T)

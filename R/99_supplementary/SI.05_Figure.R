# SCRIPT:          figure3.R
# AUTHOR:          G.B. (Giaco)
# PURPOSE:         Plot supplementary figure

rm(list = ls())

library(tidyverse)   # Data wrangling
library(patchwork)   # Plot wrangling

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

# FUNCTIONS ####
source(sprintf("%s/%s/function/iam.date.fun.plot.R", wd_oa, wd_oa_scripts))

# DATA PREPARATION ####
fam_cor <- read_csv(sprintf("%s/%s/SI.03_fam_cor_EA.csv", wd_noa, wd_noa_data)) %>%
  mutate(trait = "EA") %>%
  dplyr::select(-c(name, n)) %>%
  mutate(correlation = "Observed", title = "Resemblance between relatives") %>%
  add_row(type = "ram_s", lbound = NA, estimate = NA, ubound = NA,
          trait = "EA", correlation = "Observed",
          title = "Resemblance between relatives")

est <- read.csv(sprintf("%s/SI/SI.04_table_est_EA.csv", wd_oa)) %>%
  mutate(trait = "EA")

# Best fitting model
mod_best <- "iAM-ACE"
mod_comp <- "dAM-ACE"

# Implied spousal correlation on S (for ram_s annotation)
ram_s <- est %>%
  filter(model == mod_comp | model == mod_best , parameter == "r_Partner_S") %>%
  dplyr::select(model, estimate, lbound, ubound)

implied <- est %>%
  filter(str_starts(parameter, "r")) %>%
  mutate(
    type       = recode(parameter, !!!recode_type),
    model      = recode(model,
                        "implied (iAM-ACE)" = mod_best,
                        "implied (dAM-ACE)" = mod_comp),
    title      = "Resemblance between relatives",
    lbound     = NA_real_,
    ubound     = NA_real_
  ) %>%
  dplyr::select(type, lbound, estimate, ubound, trait, model, title) %>%
  filter(type %in% type_levels) %>%
  rename(correlation = model)

fam_cor_long <- rbind(fam_cor, implied) %>%
  mutate(
    lbound   = ifelse(type == "ram_s" & correlation == mod_comp, ram_s %>% filter(model == mod_comp) %>% pull(lbound),   lbound),
    estimate = ifelse(type == "ram_s" & correlation == mod_comp, ram_s %>% filter(model == mod_comp) %>% pull(estimate), estimate),
    ubound   = ifelse(type == "ram_s" & correlation == mod_comp, ram_s %>% filter(model == mod_comp) %>% pull(ubound),   ubound),
    lbound   = ifelse(type == "ram_s" & correlation == mod_best, ram_s %>% filter(model == mod_best) %>% pull(lbound),   lbound),
    estimate = ifelse(type == "ram_s" & correlation == mod_best, ram_s %>% filter(model == mod_best) %>% pull(estimate), estimate),
    ubound   = ifelse(type == "ram_s" & correlation == mod_best, ram_s %>% filter(model == mod_best) %>% pull(ubound),   ubound)
  )

var <- est %>%
  filter(model == mod_best,
         parameter %in% c("varP_A","varP_C","varP_E","varS_A","varS_C","varS_E")) %>%
  mutate(
    trait = case_when(
      str_starts(parameter, "varP") ~ "Focal trait",
      str_starts(parameter, "varS") ~ "Sorting factor (S)"
    ),
    parameter = case_when(
      parameter %in% c("varP_A","varS_A") ~ "h2",
      parameter %in% c("varP_C","varS_C") ~ "c2",
      parameter %in% c("varP_E","varS_E") ~ "e2"
    ),
    trait = factor(trait, levels = c("Focal trait","Sorting factor (S)"))
  )

# PLOT ####
panel_a <- make_panel_a(
  fam_cor_long  = fam_cor_long,
  implied_label = "iAM-ACE",
  implied_color = "#E59F00",
  fill_values   = c("iAM-ACE" = "#E59F00", "dAM-ACE" = "lightgray", "Observed" = "#0c72b2"),
  shape_values  = c(24, 22, 21),
  y_max = 1.1,
  trait_name = "EA"
)

panel_b <- make_panel_b(
  est          = var,
  param_levels = c("h2","c2","e2"),
  x_labels     = c(
    h2 = expression(hat(italic(h)^{"2"})),
    c2 = expression(hat(c^{"2"})),
    e2 = expression(hat(e^{"2"}))
  ),
  fill_colors  =  c("#b41f01","#008B8B","#3e59dc","#b41f01","#008B8B","#3e59dc")
)

# SAVE ####
plot <- (panel_a + panel_b) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = sprintf("%s/%s/SI.05.FIG_3.pdf", wd_oa, wd_oa_figures),
  plot     = plot,
  width    = 12.5,
  height   = 3.25,
  device   = cairo_pdf
)
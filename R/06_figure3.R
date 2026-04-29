# SCRIPT:          figure3.R
# AUTHOR:          G.B. (Giaco)
# PURPOSE:         Plot figure 3

rm(list = ls())

library(tidyverse)   # Data wrangling
library(ggokabeito)  # Colour for plotting

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
fam_cor_ht <- read_csv(sprintf("%s/%s/02.04_fam_cor_height.csv", wd_noa, wd_noa_data)) %>%
  dplyr::select(-c(name,n)) %>% 
  mutate(Trait = "Height") 
fam_cor_oe <- read_csv(sprintf("%s/%s/03.04_fam_cor_openness.csv", wd_noa, wd_noa_data)) %>%
  dplyr::select(-c(name,n)) %>% 
  mutate(Trait = "Openness") 
fam_cor_bp <- read_csv(sprintf("%s/%s/04.04_fam_cor_bpd.csv", wd_noa, wd_noa_data)) %>%
  dplyr::select(-c(name,n)) %>% 
  mutate(Trait = "Borderline") 
fam_cor_ea <- read_csv(sprintf("%s/%s/SI.03_fam_cor_EA.csv", wd_noa, wd_noa_data)) %>%
  dplyr::select(-c(name,n)) %>% 
  mutate(Trait = "Educational\nattainment") 

trait_levels <- c("Height", "Openness", "Borderline", "Educational\nattainment")

est_ea <- read.csv(sprintf("%s/SI/SI.04_table_est_EA.csv", wd_oa)) %>%
  mutate(Trait = "Educational\nattainment") %>% 
  filter(startsWith(model,"dAM"))%>%
  mutate(
    type       = recode(parameter, !!!recode_type),
    lbound     = NA_real_,
    ubound     = NA_real_,
    correlation = "Implied"
  ) %>%
  dplyr::select(type, lbound, estimate, ubound, Trait, correlation) %>%
  filter(type %in% type_levels & !type %in% c("ram", "ram_s"))

est_ht <- read.csv(sprintf("%s/SI/02.06_table_est_height.csv", wd_oa)) %>%
  mutate(Trait = "Height") %>% 
  filter(startsWith(model,"dAM"))%>%
  mutate(
    type       = recode(parameter, !!!recode_type),
    lbound     = NA_real_,
    ubound     = NA_real_,
    correlation = "Implied"
  ) %>%
  dplyr::select(type, lbound, estimate, ubound, Trait, correlation) %>%
  filter(type %in% type_levels & !type %in% c("ram", "ram_s"))

est_oe <- read.csv(sprintf("%s/SI/03.06_table_est_openness.csv", wd_oa)) %>%
  mutate(Trait = "Openness") %>% 
  filter(startsWith(model,"dAM"))%>%
  mutate(
    type       = recode(parameter, !!!recode_type),
    lbound     = NA_real_,
    ubound     = NA_real_,
    correlation = "Implied"
  ) %>%
  dplyr::select(type, lbound, estimate, ubound, Trait, correlation) %>%
  filter(type %in% type_levels & !type %in% c("ram", "ram_s"))

est_bp <- read.csv(sprintf("%s/SI/04.06_table_est_bpd.csv", wd_oa)) %>%
  mutate(Trait = "Borderline") %>% 
  filter(startsWith(model,"dAM"))%>%
  mutate(
    type       = recode(parameter, !!!recode_type),
    lbound     = NA_real_,
    ubound     = NA_real_,
    correlation = "Implied"
  ) %>%
  dplyr::select(type, lbound, estimate, ubound, Trait, correlation) %>%
  filter(type %in% type_levels & !type %in% c("ram", "ram_s"))

# save crosses values for supplementary
write_csv(est_ea %>% dplyr::select(type, estimate) %>% mutate(estimate = round(estimate,4)), sprintf("%s/SI/08_table_implied_cor_EA.csv", wd_oa))
write_csv(est_ht %>% dplyr::select(type, estimate) %>% mutate(estimate = round(estimate,4)), sprintf("%s/SI/08_table_implied_cor_height.csv", wd_oa))
write_csv(est_oe %>% dplyr::select(type, estimate) %>% mutate(estimate = round(estimate,4)), sprintf("%s/SI/08_table_implied_cor_openness.csv", wd_oa))
write_csv(est_bp %>% dplyr::select(type, estimate) %>% mutate(estimate = round(estimate,4)), sprintf("%s/SI/08_table_implied_cor_bpd.csv", wd_oa))


df_plot <- rbind(fam_cor_ht, fam_cor_oe, fam_cor_bp, fam_cor_ea) %>%
  mutate(
    Trait = factor(Trait, levels = trait_levels),
    type  = factor(type, levels = type_levels),
    correlation = "Observed"
    ) %>%
  rbind(est_ea, est_ht, est_oe, est_bp) %>%
  filter(!type == "rdz-in-law") 

# PLOT ####
 plot <- ggplot(df_plot %>% filter(Trait != "Educational\nattainment"), aes(x = type, y = estimate, fill = Trait)) +
  geom_vline(xintercept = 1.5, color = "gray", linetype = "dashed") +
  geom_hline(yintercept = 0,   color = "lightgray") +
  geom_linerange(
    aes(ymin = lbound, ymax = ubound),
    color = "gray", size = 0.75, position = position_dodge(0.5)
  ) +
  geom_point(data = df_plot %>% filter(correlation == "Observed" & Trait != "Educational\nattainment"), shape = 21, size = 2.5, position = position_dodge(0.5)) +
   scale_fill_manual(
     name = "Trait",
     values = c(unname(palette.colors(palette = "Okabe-Ito")[2:4]),"lightgray")
   ) +
  geom_point(data = df_plot %>% filter(correlation == "Implied" & Trait != "Educational\nattainment"), aes(color = Trait, shape = correlation), size = 1.5, stroke = 1, position = position_dodge(0.5)) +
  scale_y_continuous(limits = c(-.2, 1), breaks = seq(-.2, 1, 0.2)) +
  scale_shape_manual(name = "correlation", values = c(4), guide = "none") +
   scale_color_manual(
     name = "Trait",
     values = c(unname(palette.colors(palette = "Okabe-Ito")[2:4]),"lightgray")
   ) +
  scale_x_discrete(labels = function(x) str_wrap(type_labels[x], width = 8)) +
  labs(y = "Estimate", x = "Relationship type") +
  theme_bw_classic() +
  theme()

 plot_si <- ggplot(df_plot %>% filter(Trait == "Educational\nattainment"), aes(x = type, y = estimate)) +
   geom_vline(xintercept = 1.5, color = "gray", linetype = "dashed") +
   geom_hline(yintercept = 0,   color = "lightgray") +
   geom_linerange(
     aes(ymin = lbound, ymax = ubound),
     color = "gray", size = 0.75, position = position_dodge(0.5)
   ) +
   geom_point(data = df_plot %>% filter(correlation == "Observed" & Trait == "Educational\nattainment"), shape = 21, size = 2.5, position = position_dodge(0.5), fill = "lightgray") +
   geom_point(data = df_plot %>% filter(correlation == "Implied" & Trait == "Educational\nattainment"), aes(shape = correlation), size = 1.5, stroke = 1, position = position_dodge(0.5),  color = "lightgray") +
   scale_y_continuous(limits = c(-.2, 1), breaks = seq(-.2, 1, 0.2)) +
   scale_shape_manual(name = "correlation", values = c(4), guide = "none") +
   scale_x_discrete(labels = function(x) str_wrap(type_labels[x], width = 8)) +
   labs(y = "Estimate", x = "Relationship type") +
   theme_bw_classic() +
   theme()
 
# SAVE ####
ggsave(
  filename = sprintf("%s/%s/FIG_3.pdf", wd_oa, wd_oa_figures),
  plot     = plot,
  width    = 8,
  height   = 3,
  device   = cairo_pdf
)

ggsave(
  filename = sprintf("%s/%s/06.FIG_3.S.pdf", wd_oa, wd_oa_figures),
  plot     = plot_si,
  width    = 8,
  height   = 3,
  device   = cairo_pdf
)

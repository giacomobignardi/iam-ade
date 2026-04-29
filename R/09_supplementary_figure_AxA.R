# SCRIPT:          figure3.R
# AUTHOR:          G.B. (Giaco)
# PURPOSE:         Plot figure 3

rm(list = ls())

library(tidyverse)   # Data wrangling
library(ggokabeito)  # Colour for plotting
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

# Define trait order for consistent plotting
trait_levels <- c("Educational\nattainment", "Height", "Openness", "Borderline")

# DATA PREPARATION ####
# Helper function to load and prepare trait estimates
load_trait_estimates <- function(file_path, trait_name) {
  data <- read.csv(file_path) %>%
    mutate(Trait = trait_name)
 
  # data <- data %>% filter(axa %in% c(2,4,6,8,10))
  # 
  data %>%
    dplyr::select(axa, rNA, parameter, lbound, estimate, ubound, Trait, model)
}

# Load ACE model estimates (trait-specific models)
est_ht <- load_trait_estimates(
  sprintf("%s/SI/02.06_table_est_nonadditive_height.csv", wd_oa),
  "Height"
)

est_oe <- load_trait_estimates(
  sprintf("%s/SI/03.06_table_est_nonadditive_openness.csv", wd_oa),
  "Openness"
)

est_bp <- load_trait_estimates(
  sprintf("%s/SI/04.06_table_est_nonadditive_bpd.csv", wd_oa),
  "Borderline"
)

# Parameter recoding: convert internal names to publication labels
parameter_recode <- c(
  "A" = "A",  # Additive variance
  "N" = "I",  # Dominance variance
  "C" = "C",  # Common environmental variance
  "E" = "E"   # Unique environmental variance
)

# Combine ACE model estimates and recode parameters
est <- rbind(est_ht, est_oe, est_bp) %>%
  filter(parameter %in% names(parameter_recode)) %>%
  mutate(parameter = recode(parameter, !!!parameter_recode),
         parameter = factor(parameter, levels = c("A", "I", "E")),
         Trait = factor(Trait, levels = trait_levels)
  )

est %>% filter(parameter %in% c("A", "I") & Trait == "Height")  # 4 and 5
est %>% filter(parameter %in% c("A", "I") & Trait == "Openness") # 4 and 6
est %>% filter(parameter %in% c("A", "I") & Trait == "Borderline")  # 4 and 5

est %>% filter(parameter %in% c("A", "I") & axa %in% c(2,10))

# PLOT ####
panel_a <- ggplot(
  est %>% filter(estimate != 0 ),  # Exclude zero estimates
  aes(x = as.numeric(axa), y = estimate, color = parameter, fill = parameter)
) +
  # Line trends for each parameter
  geom_line(
    linewidth = 1,
    aes(linetype = parameter)
  ) +
  # Points at data locations
  geom_point(
    size = 3,
    shape = 21,
    color = "black",
    stroke = 0.8
  ) +
  # Confidence intervals as ribbons
  geom_ribbon(
    aes(ymin = lbound, ymax = ubound, fill = parameter),
    alpha = 0.2,
    color = NA
  ) +
  # Y-axis formatting
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0),
    name = "(Aprox.) Proportion of variance\n under epistasis"
  ) +
  # X-axis: numeric for line trends
  scale_x_continuous(
    breaks = unique(as.numeric(est$axa)),
    name = "AxA interaction order"
  ) +
  # Color palette
  scale_color_manual(
    values = c(
      A = "#b41f01",  # Red
      I = "#d87a38",  # Orange
      C = "#008B8B",  # Teal
      E = "#3e59dc"   # Blue
    ),
    name = "Parameter"
  ) +
  # Fill for ribbons
  scale_fill_manual(
    values = c(
      A = "#b41f01",
      I = "#d87a38",
      C = "#008B8B",
      E = "#3e59dc"
    ),
    guide = "none"
  ) +
  # Line types for distinction
  scale_linetype_manual(
    values = c(A = "solid", I = "dashed", C = "dotted", E = "twodash"),
    guide = "none"
  ) +
  # Facets
  facet_wrap(~Trait, nrow = 1, scales = "free_x") +
  # Theme
  theme_bw_classic() +
  theme(legend.position = "bottom")

# Display plot
panel_a

# SAVE ####
ggsave(
  filename = sprintf("%s/%s/09.FIG_S7.pdf", wd_oa, wd_oa_figures),
  plot     = panel_a,
  width    = 8,
  height   = 4,
  device   = cairo_pdf
)
est %>% filter(axa %in% c(2,10) )

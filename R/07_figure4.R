# SCRIPT:          figure3.R
# AUTHOR:          G.B. (Giaco)
# PURPOSE:         Plot figure 3

rm(list = ls())

library(tidyverse)   # Data wrangling

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
trait_levels <- c("Height", "Openness", "Borderline","Educational\nattainment")

# DATA PREPARATION ####
# Helper function to load and prepare trait estimates
load_trait_estimates <- function(file_path, trait_name, model_filter = NULL) {
  data <- read.csv(file_path) %>%
    mutate(Trait = trait_name)
  
  if (!is.null(model_filter)) {
    data <- data %>% filter(model == model_filter)
  }
  
  data %>%
    dplyr::select(parameter, lbound, estimate, ubound, Trait, model)
}

# Load ACE model estimates (trait-specific models)
est_ea <- load_trait_estimates(
  sprintf("%s/SI/SI.04_table_est_EA.csv", wd_oa),
  "Educational\nattainment",
  "iAM-ACE"
)

est_ht <- load_trait_estimates(
  sprintf("%s/SI/02.06_table_est_height.csv", wd_oa),
  "Height",
  "xdAM-ADE"
)

est_oe <- load_trait_estimates(
  sprintf("%s/SI/03.06_table_est_openness.csv", wd_oa),
  "Openness",
  "dAM-ADE"
)

est_bp <- load_trait_estimates(
  sprintf("%s/SI/04.06_table_est_bpd.csv", wd_oa),
  "Borderline",
  "xdAM-ADE"
)

# Load CTD (Cross-Twin Differences) model estimates (all models included)
est_ea_ctd <- load_trait_estimates(
  sprintf("%s/SI/SI.04_tableCTD_est_EA.csv", wd_oa),
  "Educational\nattainment"
)

est_ht_ctd <- load_trait_estimates(
  sprintf("%s/SI/02.06_tableCTD_est_height.csv", wd_oa),
  "Height"
)

est_oe_ctd <- load_trait_estimates(
  sprintf("%s/SI/03.06_tableCTD_est_openness.csv", wd_oa),
  "Openness"
)

est_bp_ctd <- load_trait_estimates(
  sprintf("%s/SI/04.06_tableCTD_est_bpd.csv", wd_oa),
  "Borderline"
)

# Parameter recoding: convert internal names to publication labels
parameter_recode <- c(
  "varP_A" = "A",  # Additive variance
  "varP_D" = "D",  # Dominance variance
  "varP_C" = "C",  # Common environmental variance
  "varP_E" = "E"   # Unique environmental variance
)

# Combine ACE model estimates and recode parameters
est <- rbind(est_ea, est_ht, est_oe, est_bp) %>%
  filter(parameter %in% names(parameter_recode)) %>%
  mutate(parameter = recode(parameter, !!!parameter_recode))

# Combine CTD model estimates (exclude V parameter)
est_ctd <- rbind(est_ea_ctd, est_ht_ctd, est_oe_ctd, est_bp_ctd) %>%
  filter(parameter != "V")

# Merge all estimates and set factor levels for consistent ordering
est <- rbind(est, est_ctd) %>%
  mutate(
    parameter = factor(parameter, levels = c("A", "D", "C", "E")),
    Trait = factor(Trait, levels = trait_levels)
  )

# PLOT ####

shading <- data.frame(
  Trait = "Educational\nattainment",  # must match your facet label exactly
  xmin = -Inf, xmax = Inf,
  ymin = -Inf, ymax = Inf
)


panel_a <- ggplot(
  est %>% filter(estimate != 0 & Trait != "Educational\nattainment"),  # Exclude zero estimates (empty bars)
  aes(x = model, y = estimate, fill = parameter)
)  +
  # Bars with consistent width across facets
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8, preserve = "single"),
    color = "black",
    width = 0.7
  ) +
  # Confidence intervals
  geom_errorbar(
    aes(ymin = lbound, ymax = ubound),
    position = position_dodge(width = 0.8, preserve = "single"),
    width = 0.25,
    linewidth = 0.5,
    color = "black"
  ) +
  # Estimate labels on bars
  geom_text(
    aes(label = round(estimate, 2), y = ubound),
    position = position_dodge(width = 0.8, preserve = "single"),
    vjust = -0.8,
    size = 3
  ) +
  # Y-axis formatting
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0),
    name = "Proportion of variance"
  ) +
  # X-axis with mathematical notation for parameters
  scale_x_discrete(
    labels = c(
      A = expression(hat(italic(h)^{"2"})),     # Heritability
      D = expression(hat(italic(delta)^{"2"})), # Dominance variance
      C = expression(hat(italic(c)^{"2"})),     # Common environment
      E = expression(hat(italic(e)^{"2"}))      # Unique environment
    ),
    drop = TRUE,  # Drop unused model levels per trait
    name = "Most parsimonious model"
  ) +
  # Color palette for variance components
  scale_fill_manual(
    values = c(
      A = "#b41f01",  # Red
      D = "#faa53c",  # Orange
      C = "#008B8B",  # Teal
      E = "#3e59dc"   # Blue
    )
  ) +
  # Separate panels for each trait with independent x-axes
  facet_wrap(~Trait, nrow = 1, scales = "free_x") +
  # Theme and styling
  theme_bw_classic()

# Display plot
panel_a


ggsave(
  filename = sprintf("%s/%s/FIG_4.pdf", wd_oa, wd_oa_figures),
  plot     = panel_a,
  width    = 8,
  height   = 3,
  device   = cairo_pdf
)


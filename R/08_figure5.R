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

# DATA ####
trait_levels <- c("Educational\nattainment", "Height", "Openness", "Borderline")

# load estimates
est_ht <- read.csv(sprintf("%s/SI/02.06_table_est_height.csv", wd_oa)) %>%
  mutate(Trait = "Height",
         Model = "Bt+Wt variance") %>% 
  filter(model == "xdAM-ADE") %>% 
  dplyr::select(parameter, lbound, estimate, ubound, Trait, Model) 
est_oe <- read.csv(sprintf("%s/SI/03.06_table_est_openness.csv", wd_oa)) %>%
  mutate(Trait = "Openness",
         Model = "Bt+Wt variance") %>% 
  filter(model == "dAM-ADE") %>% 
  dplyr::select(parameter, lbound, estimate, ubound, Trait, Model) 
est_bp <- read.csv(sprintf("%s/SI/04.06_table_est_bpd.csv", wd_oa)) %>%
  mutate(Trait = "Borderline",
         Model = "Bt+Wt variance") %>% 
  filter(model == "xdAM-ADE") %>% 
  dplyr::select(parameter, lbound, estimate, ubound, Trait, Model) 

# load estimates adjusted for within-individual variance
est_ht_err <- read.csv(sprintf("%s/SI/02.06_table_est_error_height.csv", wd_oa)) %>%
  mutate(Trait = "Height",
         Model = "Bt variance") %>% 
  dplyr::select(parameter, lbound, estimate, ubound, Trait, Model) 
est_oe_err <- read.csv(sprintf("%s/SI/03.06_table_est_error_openness.csv", wd_oa)) %>%
  mutate(Trait = "Openness",
         Model = "Bt variance") %>% 
  dplyr::select(parameter, lbound, estimate, ubound, Trait, Model) 
est_bp_err <- read.csv(sprintf("%s/SI/04.06_table_est_error_bpd.csv", wd_oa)) %>%
  mutate(Trait = "Borderline",
         Model = "Bt variance") %>% 
  dplyr::select(parameter, lbound, estimate, ubound, Trait, Model) 

# load ICC
load(sprintf("%s/%s/02.01_rYY_height.Rdata", wd_noa, wd_noa_data))
ryy_ht <- rYY$ICC
load(sprintf("%s/%s/03.01_rYY_openness.Rdata", wd_noa, wd_noa_data))
ryy_oe <- rYY$ICC
load(sprintf("%s/%s/04.01_rYY_bpd.Rdata", wd_noa, wd_noa_data))
ryy_bp <- rYY$ICC

recode_type <- c(
  "varP_A"   = "A",
  "varP_D"   = "D",
  "r_Partner_P" = "r_mate",
  "rA_Partner"  = "rA_mate",
  "rA_FS"       = "rA_FS"
  
)

recode_type_adj <- c(
  "varP_bw_A"   = "A",
  "varP_bw_D"   = "D",
  "r_Partner_S" = "r_mate (bw)",
  "rA_Partner"  = "rA_mate (bw)",
  "rA_FS"       = "rA_FS"
  
)

est_unadj <- rbind(est_ht, est_oe, est_bp) %>% 
  filter(parameter %in% c("varP_A","varP_D", "r_Partner_P", "rA_Partner" , "rA_FS")) %>% 
  mutate(parameter = recode(parameter, !!!recode_type))

est_adj <- rbind(est_ht_err, est_oe_err, est_bp_err)%>% 
  filter(parameter %in% c("varP_bw_A","varP_bw_D", "r_Partner_S", "rA_Partner" , "rA_FS")) %>% 
  mutate(parameter  = recode(parameter, !!!recode_type_adj))

est <- rbind(est_unadj, est_adj) %>% 
  mutate(  
    parameter = factor(parameter, levels = c("A", "D", "r_mate", "rA_mate", "r_mate (bw)", "rA_mate (bw)", "rA_FS")),
    Trait     = factor(Trait    , levels = c("Height", "Openness", "Borderline")),
    Model     = factor(Model    , levels = c("Bt+Wt variance", "Bt variance"))
  )

# Create summary for totals
est_broad <- est %>% 
  filter(parameter %in% c("A","D")) %>%
  group_by(Trait, Model) %>%
  summarise(broad = sum(estimate))

# PLOT ####
# plot broad sense heritabilits
panel_a <- ggplot(est %>% filter(parameter %in% c("A","D")), aes(x = Trait, y = estimate)) +
  geom_bar(aes(fill = parameter), stat = "identity", position = "stack", width = .7, color = "black") +
  geom_text(
    aes(fill = parameter, label = round(estimate,2)),
    position = position_stack(vjust = 0.5),  # centers text within each segment
    color = "white",
    size = 3
  ) +
  geom_text(
    data = est_broad,
    aes(x = Trait, y = broad, label = sprintf("%.2f", broad)),
    vjust = -0.5,  
    color = "black",
    size = 3
  ) +
  facet_wrap(~Model, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1, 0.2)) +
  scale_fill_manual(
    values = c("A" = "#b41f01", "D" = "#faa53c")
  ) +
  labs(
    x = "Trait",
    y = "Proportion of variance",
    fill = "Parameter"
  ) +
  theme_bw_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

  panel_b <- ggplot(est %>% filter(!parameter %in% c("A","D")), aes(x = parameter, y = estimate, fill = Trait)) +
  geom_linerange(
    aes(ymin = lbound, ymax = ubound),
    color = "gray", size = 0.75, position = position_dodge(0.5)
  ) +
     geom_errorbar(
      data = est %>% filter(parameter == "rA_FS") %>% distinct(Model) %>% mutate(parameter = "rA_FS"),
      aes(x = parameter, ymin = 0.5, ymax = 0.5),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed"
    ) +
    geom_errorbar(
      data = est %>% filter(parameter %in% c("rA_mate","rA_mate (bw)")) %>% distinct(Model) %>% mutate(parameter = c("rA_mate","rA_mate (bw)")),
      aes(x = parameter, ymin = 0, ymax = 0),
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed"
    ) +
  facet_wrap(~Model, nrow = 1, scales = "free_x") +
  geom_point(shape = 21, size = 2.5, position = position_dodge(0.5)) +
  scale_fill_okabe_ito(name = "Trait") +
  scale_y_continuous(limits = c(0, .75), breaks = seq(0,.8, 0.2)) +
  labs(y = "Estimate", x = "Relationship type") +
  scale_x_discrete(labels =
                     c(
    "r_mate" = "spouse\n(focal)",
    "rA_mate" = "spouse\n(focal)\ngenetic",
    "r_mate (bw)" = "spouse\n(sorting)",
    "rA_mate (bw)" = "spouse\n(sorting)\ngenetic",
    "rA_FS" ="sibling|DZ\ngenetic")
    ) +
  theme_bw_classic() +
  theme()

plot <- panel_a + panel_b + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect", widths = c(1.5, 2))

ggsave(
  filename = sprintf("%s/%s/FIG_5.pdf", wd_oa, wd_oa_figures),
  plot     = plot,
  width    = 10,
  height   = 3,
  device   = cairo_pdf
)


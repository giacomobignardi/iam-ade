# SCRIPT:  simulations.R
# AUTHOR:  G.B. (Giaco)
# PURPOSE: Simulate twin and extended family data under various assortative
#          mating (AM) scenarios. Computes expected phenotypic correlations
#          across family relationships (PART 1), and generates raw simulated
#          datasets for model fitting via CTD/ETFD (PART 2).

rm(list = ls())

library(tidyverse)  # Data wrangling and plotting
library(OpenMx)     # Structural equation modelling
library(MASS)       # Multivariate data simulation
library(patchwork)  # Panle plotting wrangling

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


# FUNCTIONS ####
# Simulate standardised (expected) family correlations under a given AM model
source(sprintf("%s/%s/function/sim.estd.fun.R",      wd_oa, wd_oa_scripts))

# Classical Twin Design (CTD) model: ACE / ADE decomposition
source(sprintf("%s/%s/function/acde.fun.R",          wd_oa, wd_oa_scripts))

# Extended Twin Family Design (ETFD) / IAM-DATE model
source(sprintf("%s/%s/function/iam.date.fun.R",  wd_oa, wd_oa_scripts))

# Plotting helpers for IAM-DATE figures
source(sprintf("%s/%s/function/iam.date.fun.plot.R", wd_oa, wd_oa_scripts))

# PART 1: EXPECTED CORRELATIONS ACROSS FAMILY RELATIONSHIPS ####
# Sweeps over:
#- Assortment type
#- Magnitude of assortment on the sorting factor (mu)
#- Reliability of the measure (ryy)
#- Preferential mating (xi)
#- Stratification (q2).

# Strength of assortative mating (co-path coefficient mu)
mu.vals  <- seq(0, .9, .1)

# Variance component proportions (squared path coefficients; must sum to <= 1)
a2.val <- .30  # Additive genetic variance
d2.val <- .15  # Dominance genetic variance
c2.val <- .00  # Shared (common) environmental variance
i2.val <- .00  # Non-additive genetic variance beyond dominance (e.g. epistasis)
n2.val <- d2.val + i2.val  # Total non-additive genetic variance
t2.val <- .00  # Twin-shared social homogamy variance

# Measurement reliability (phenotypic, between-individual) capped at .5 (<.5 deemed to be poor)
ryy.vals <- seq(.5, 1, .1)

# Prefential mating as cross-trait assortment × sibling shared mate preferences (interaction term) capped at .5 (arbitrary choice)
xi.vals  <- seq(0, .5, .1)

# Types of assortative mating modelled
am.type <- c(
  "dAM",               # Direct (phenotypic) assortative mating
  "narrow-sense hom.", # Homogamy based on narrow-sense h2
  "broad-sense hom.",  # Homogamy based on broad-sense H2
  "idio. hom.",        # Idiosyncratic (individual-specific) homogamy
  "mixed hom."         # Mixture
)

# Display labels for family-relationship correlation components
component_labels <- c(
  "rmz"    = "MZ",
  "rsib"   = "sibling|DZ",
  "rmzil"  = "MZ in-law",
  "rdzil"  = "sibling|DZ in-law",
  "rmzcil" = "MZ co-in-law",
  "rdzcil" = "DZ co-in-law",
  "ram"    = "spouse (focal)"
)

# NESTED SIMULATION ####
# Iterates over: reliability × xi × stratification × AM type × mu
# Results are stored in a 5-level nested named list.
sim_results <- setNames(vector("list", length(ryy.vals)), as.character(ryy.vals))

for (ryy.val in ryy.vals) {
  
  sim_results[[as.character(ryy.val)]] <-
    setNames(vector("list", length(xi.vals)), as.character(xi.vals))
  
  for (xi.val in xi.vals) {
    
    # Stratification grid (proportion of variance due to population structure)
    q2.vals <- seq(.0, .5, .1)
    
    sim_results[[as.character(ryy.val)]][[as.character(xi.val)]] <-
      setNames(vector("list", length(q2.vals)), as.character(q2.vals))
    
    for (q2.val in q2.vals) {
      
      sim_results[[as.character(ryy.val)]][[as.character(xi.val)]][[as.character(q2.val)]] <-
        setNames(vector("list", length(am.type)), am.type)
      
      for (type in am.type) {
        
        sim_results[[as.character(ryy.val)]][[as.character(xi.val)]][[as.character(q2.val)]][[type]] <-
          setNames(vector("list", length(mu.vals)), as.character(mu.vals))
        
        for (mu.val in mu.vals) {
          
          # Call sim.estd(); wrap in tryCatch so a single failure does not
          # abort the full sweep — failed cells are stored as NULL.
          sim_results[[as.character(ryy.val)]][[as.character(xi.val)]][[as.character(q2.val)]][[type]][[as.character(mu.val)]] <-
            tryCatch(
              {
                sim.estd(
                  a2.val  = a2.val,
                  d2.val  = d2.val,
                  i2.val  = i2.val,
                  c2.val  = c2.val,
                  t2.val  = t2.val,
                  mu.val  = mu.val,
                  ryy     = ryy.val,
                  xi      = xi.val,
                  q2      = q2.val,
                  Nmz     = 5000,    # MZ twin pairs (sample stats true - sample invariant)
                  Ndz     = 5000,    # DZ twin pairs (sample stats true - sample invariant)
                  am.type = type
                )
              },
              error = function(e) {
                message(sprintf(
                  "Skipping: ryy = %s, xi = %s, q2 = %s, type = %s, mu = %s | Error: %s",
                  ryy.val, xi.val, q2.val, type, mu.val, e$message
                ))
                return(NULL)
              }
            )
        }
      }
    }
  }
}

# EXTRACT RESULTS (COR) ####
# Identify the first non-NULL result to use as a column-name template
first_result <- NULL
for (ryy.val in ryy.vals) {
  for (xi.val in xi.vals) {
    for (q2.val in q2.vals) {
      for (type in am.type) {
        for (mu.val in mu.vals) {
          result <- sim_results[[as.character(ryy.val)]][[as.character(xi.val)]][[as.character(q2.val)]][[type]][[as.character(mu.val)]]
          if (!is.null(result)) { first_result <- result; break }
        }
        if (!is.null(first_result)) break
      }
      if (!is.null(first_result)) break
    }
    if (!is.null(first_result)) break
  }
  if (!is.null(first_result)) break
}

# Initialise output data frame with correct column structure (zero rows)
est_am <- as.data.frame(first_result)[0, ]

# Flatten the nested list, tagging each row with its parameter combination
for (ryy.val in ryy.vals) {
  for (xi.val in xi.vals) {
    for (q2.val in q2.vals) {
      for (type in am.type) {
        for (mu.val in mu.vals) {
          
          result <- sim_results[[as.character(ryy.val)]][[as.character(xi.val)]][[as.character(q2.val)]][[type]][[as.character(mu.val)]]$value
          
          if (is.null(result)) next  # Skip failed / missing cells
          
          df            <- as.data.frame(result)
          df$ryy        <- ryy.val
          df$xi         <- xi.val
          df$q2         <- q2.val
          df$assortment <- type
          df$mu         <- mu.val
          
          est_am <- rbind(est_am, df)
        }
      }
    }
  }
}

# Reshape to long format for ggplot
est_am_long <- est_am %>%
  as.data.frame() %>%
  pivot_longer(
    cols      = -c(assortment, mu, q2, xi, ryy),
    names_to  = "component",
    values_to = "est"
  ) %>%
  mutate(across(everything(), unlist)) %>%
  # Apply human-readable labels defined in component_labels
  mutate(component = recode(component, !!!component_labels))


# FIGURE 2 — Panel A ####
# Main result: no stratification (q2 = 0), perfect reliability (ryy = 1),
# no preferential mating (xi = 0).
# Colour encodes partner correlation (mu).
plot_a <- est_am_long %>%
  filter(q2 == 0, ryy == 1, xi == 0) %>%
  filter(!component %in% c("Vp", "Vs", "rsibil")) %>%
  mutate(
    component  = factor(component,  levels = rev(component_labels)),
    assortment = factor(assortment, levels = c(
      "dAM", "narrow-sense hom.", "broad-sense hom.",
      "social hom.", "idio. hom.", "mixed hom."
    ))
  ) %>%
  ggplot(aes(assortment, est, fill = mu)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_gradientn(colours = c(
    "#120246", 
    "#050353", 
    "#06175F", 
    "#08306B", 
    "#285C80",
    "#478495", 
    "#68A7A9", 
    "#89BDB5", 
    "#AAD0C4", 
    "#CCE3D7"
  )) +
  facet_grid(cols = vars(component), scales = "free") +
  labs(
    x    = expression(paste("Causes of partner correlation")),
    y    = expression(paste("Phenotypic correlation")),
    fill = expression(mu)
  ) +
  ylim(c(-.01, 1)) +
  theme_bw_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1)
  )

## FIGURE SI — S1 #### 
# Variation over measurement reliability (ryy)
# Fixed: q2 = 0, mu = 0.5, xi = 0. Colour encodes ryy.
plot_a.s1 <- est_am_long %>%
  filter(q2 == .0, mu == .5, xi == 0) %>%
  filter(!component %in% c("Vp", "Vs", "rsibil")) %>%
  mutate(
    component  = factor(component,  levels = rev(component_labels)),
    assortment = factor(assortment, levels = c(
      "dAM", "narrow-sense hom.", "broad-sense hom.",
      "social hom.", "idio. hom.", "mixed hom."
    ))
  ) %>%
  ggplot(aes(assortment, est, fill = ryy)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_gradientn(colours = c(
    "#160830", 
    "#2A1050", 
    "#3D1A70", 
    "#562D8C", 
    "#7248A6",
    "#9168BE", 
    "#AD8FD0", 
    "#C8B0DF", 
    "#DDD0EC", 
    "#F0E8F8"
  )) +
  facet_grid(cols = vars(component), scales = "free") +
  labs(
    x    = expression(paste("Causes of partner correlation")),
    y    = expression(paste("Phenotypic correlation")),
    fill = expression(1-italic(e)[w]^2)
  ) +
  ylim(c(-.01, 1)) +
  theme_bw_classic() +
  theme(
    strip.background = element_blank(),
    
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1)
  )

## FIGURE SI — S2 ####
## Variation over preferential mating (xi)
# Fixed: q2 = 0, mu = 0.5, ryy = 1. Colour encodes xi.
plot_a.s2 <- est_am_long %>%
  filter(q2 == .0, mu == .5, ryy == 1) %>%
  filter(!component %in% c("Vp", "Vs", "rsibil")) %>%
  mutate(
    component  = factor(component,  levels = rev(component_labels)),
    assortment = factor(assortment, levels = c(
      "dAM", "narrow-sense hom.", "broad-sense hom.",
      "social hom.", "idio. hom.", "mixed hom."
    ))
  ) %>%
  ggplot(aes(assortment, est, fill = xi)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_gradientn(colours = c(
    "#021A0E", 
    "#053320", 
    "#0A4F33", 
    "#116647", 
    "#278056",
    "#469A6A", 
    "#6BB183", 
    "#96C8A0", 
    "#BDDCBC", 
    "#E0F0DC"
  )) +
  facet_grid(cols = vars(component), scales = "free") +
  labs(
    x    = expression(paste("Causes of partner correlation")),
    y    = expression(paste("Phenotypic correlation")),
    fill = expression(italic(mu)^2*italic(r)[aes])
  ) +
  ylim(c(-.01, 1)) +
  theme_bw_classic() +
  theme(
    strip.background = element_blank(),
    
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1)
  )

## FIGURE SI — S3 #### 
# Variation over population stratification (q2)
# Fixed: xi = 0, mu = 0.5, ryy = 1. Colour encodes q2.
plot_a.s3 <- est_am_long %>%
  filter(xi == .0, mu == .5, ryy == 1) %>%
  filter(!component %in% c("Vp", "Vs", "rsibil")) %>%
  mutate(
    component  = factor(component,  levels = rev(component_labels)),
    assortment = factor(assortment, levels = c(
      "dAM", "narrow-sense hom.", "broad-sense hom.",
      "social hom.", "idio. hom.", "mixed hom."
    ))
  ) %>%
  ggplot(aes(assortment, est, fill = q2)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_gradientn(colours = c(
    "#00204d", 
    "#2c4068", 
    "#898572", 
    "#c6b76f", 
    "#e9d36c"
  )) +
  facet_grid(cols = vars(component), scales = "free") +
  labs(
    x    = expression(paste("Causes of partner correlation")),
    y    = expression(paste("Phenotypic correlation")),
    fill = expression(italic(q)^2)
  ) +
  ylim(c(-.01, 1)) +
  theme_bw_classic() +
  theme(
    strip.background = element_blank(),
    
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# SAVE PART 1 ####
ggsave(
  filename = sprintf("%s/%s/00_FIG_2A.pdf", wd_oa, wd_oa_figures),
  plot     = plot_a,
  width    = 9.5,
  height   = 3,
  device   = cairo_pdf
)

ggsave(
  filename = sprintf("%s/%s/00_FIG_2A.S1.pdf", wd_oa, wd_oa_figures),
  plot     = plot_a.s1,
  width    = 9.5,
  height   = 3,
  device   = cairo_pdf
)

ggsave(
  filename = sprintf("%s/%s/00_FIG_2A.S2.pdf", wd_oa, wd_oa_figures),
  plot     = plot_a.s2,
  width    = 9.5,
  height   = 3,
  device   = cairo_pdf
)

ggsave(
  filename = sprintf("%s/%s/00_FIG_2A.S3.pdf", wd_oa, wd_oa_figures),
  plot     = plot_a.s3,
  width    = 9.5,
  height   = 3,
  device   = cairo_pdf
)


# PART 2: RAW DATA SIMULATION FOR CTD / ESTD MODEL FITTING ####
# Sweeps over additive genetic variance (a.vals) and AM type × mu
# Large N (5000 pairs, 10,000 fam) for pop. estimates.

# Parameter grid for additive heritability sweep
# Lower bound: at least twice the non-additive + social variance
# Upper bound: leaves room for residual (E) variance
a.vals <- seq(
  2 * (n2.val + t2.val),
  1  - (.1 + n2.val + t2.val),
  .05
)

# Initialise 3-level nested list: a2 × AM type × mu
sim_dat <- setNames(vector("list", length(a.vals)), as.character(a.vals))

for (a.val in a.vals) {
  
  sim_dat[[as.character(a.val)]] <-
    setNames(vector("list", length(am.type)), am.type)
  
  for (type in am.type) {
    
    sim_dat[[as.character(a.val)]][[type]] <-
      setNames(vector("list", length(mu.vals)), as.character(mu.vals))
    
    for (mu.val in mu.vals) {
      
      # Simulate raw individual-level data (raw = TRUE) with large N
      sim_dat[[as.character(a.val)]][[type]][[as.character(mu.val)]] <-
        tryCatch(
          {
            sim.estd(
              a2.val  = a.val,
              d2.val  = d2.val,
              i2.val  = i2.val,
              c2.val  = c2.val,
              t2.val  = t2.val,
              mu.val  = mu.val,
              Nmz     = 5000,   # MZ twin pairs
              Ndz     = 5000,   # DZ twin pairs
              am.type = type,
              raw     = TRUE    # Return individual-level data rather than summary stats
            )
          },
          error = function(e) {
            message(sprintf(
              "Skipping: a2 = %s, type = %s, mu = %s | Error: %s",
              a.val, type, mu.val, e$message
            ))
            return(NULL)
          }
        )
    }
  }
}


# Model fitting: CTD (ACDE) ####
# Mean structure labels — twin-only CTD uses 2 phenotypes per family;
mean.MZ     <- rep("mean", 2)
mean.DZ     <- rep("mean", 2)
sel_vars_tw <- c("tw1", "tw2") # CTD phenotype columns  

est_ctd_all <- NULL

for (a.val in a.vals) {
  for (type in am.type) {       
    for (mu.val in mu.vals) {
      
      # Retrieve simulated dataset for this parameter combination
      df <- sim_dat[[as.character(a.val)]][[type]][[as.character(mu.val)]]$dat
      
      if (is.null(df)) next       # Skip if simulation failed for this cell
      
      # Split by zygosity (1 = MZ, 2 = DZ)
      mz <- df %>% filter(zyg == 1)
      dz <- df %>% filter(zyg == 2)
      
      # Choose ACE vs ADE specification based on the Falconer criterion:
      # if 2*r_DZ - r_MZ > 0, shared environment (C) is implicated; otherwise D.
      acde_mod <- if ((2 * cor(dz$tw1, dz$tw2) - cor(mz$tw1, mz$tw2)) > 0) {
        acde.fun(mz, dz, sel_vars_tw,
                 mean.MZ = mean.MZ, mean.DZ = mean.DZ, cp.free = TRUE)  # ACE
      } else {
        acde.fun(mz, dz, sel_vars_tw,
                 mean.MZ = mean.MZ, mean.DZ = mean.MZ, dp.free = TRUE)  # ADE
      }
      
      # Fit the CTD model and extract standardised variance components
      est_ctd_i <- tryCatch(
        {
          mxRun(acde_mod, intervals = FALSE)$US$result[
            , c("varP_A", "varP_C", "varP_D", "varP_E")
          ] %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            rename(component = rowname, est = ".") %>%
            mutate(
              model      = "ACDE",
              a2         = a.val,   # Tag with generating parameter values
              assortment = type,    # so rows remain identifiable after rbind
              mu         = mu.val
            )
        },
        error = function(e) {
          message(sprintf(
            "Model fit failed: a2 = %s, type = %s, mu = %s | Error: %s",
            a.val, type, mu.val, e$message
          ))
          return(NULL)
        }
      )
      
      # Accumulate this iteration's estimates into the growing results frame
      est_ctd_all <- rbind(est_ctd_all, est_ctd_i)
      
    }
  }
}

# bias in genetic effects
est_am_bias <- est_ctd_all %>%
  # select standardized components
  filter(grepl("var", component)) %>%
  mutate(
    d2 = .15, # we simulated dominance at .15
    component = case_match(
      component,
      "varP_E" ~ "E",
      "varP_T" ~ "T",
      "varP_D" ~ "D",
      "varP_C" ~ "C",
      "varP_E" ~ "E",
      "varP_D" ~ "D",
      "varP_C" ~ "C",
      "varP_N" ~ "N",
      "varP_N" ~ "N",
      "varP_A" ~ "A",
      "varP_A" ~ "A"
    )
  ) %>%
  dplyr::filter(component %in% c("A", "D")) %>%
  pivot_wider(names_from = c("model", "component"), values_from = "est")

## Figure 2c ####
# bias in h2twin estimates
plot_b <- est_am_bias %>%
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  ggplot(aes(a2, ACDE_A, color = mu, group = mu)) +
  geom_point() +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = .30, linetype = "dashed") +
  facet_grid(cols = vars(assortment), scales = "free") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed", color = "#B41F01") +
  theme_classic(base_size = 10) +
  ylim(c(.2, 1)) +
  scale_x_continuous(breaks = seq(.30, 1, .15), limits = c(.30, 0.85)) +
  scale_color_gradientn(colors = c(
    "#120246",
    "#050353",
    "#06175F",
    "#08306B",
    "#285C80",
    "#478495",
    "#68A7A9",
    "#89BDB5",
    "#AAD0C4",
    "#CCE3D7"
  )) +
  labs(
    x = expression(italic(h)^2),
    y = expression(hat(italic(h))[twin]^2),
    color = expression(paste("Partner correlation")
                       )
  ) +
  theme_bw_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_blank()
  )

## Figure 2d ####
# bias in d2twin estimates
plot_c <- est_am_bias %>%
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  ggplot(aes(mu, ACDE_D, color = a2, group = a2)) +
  geom_hline(yintercept = .15, linetype = "dashed", color = "#FAA53C", linewidth = 1) +
  geom_point() +
  ylim(c(-0.01, 0.3)) +
  geom_line(linewidth = 1) +
  facet_grid(cols = vars(assortment), scales = "free") +
  theme_classic(base_size = 10) +
  scale_color_viridis_c(option = "plasma", breaks = seq(.3, 1, .15)) +
  scale_x_continuous(limits = c(0, .95), breaks = seq(0, 1, .20)) +
  labs(
    y = expression(hat(italic(delta))[twin]^2),
    x = expression(italic(mu)),
    color = expression(italic(h)^2)
  ) +
  theme_bw_classic(base_size = 10) +
  theme(
    strip.background = element_blank()
  )

# Initialise 3-level nested list: a2 × AM type × mu
sim_dat_pop <- setNames(vector("list", length(a.vals)), as.character(a.vals))

for (a.val in a.vals) {
  
  sim_dat_pop[[as.character(a.val)]] <-
    setNames(vector("list", length(am.type)), am.type)
  
  for (type in am.type) {
    
    sim_dat_pop[[as.character(a.val)]][[type]] <-
      setNames(vector("list", length(mu.vals)), as.character(mu.vals))
    
    for (mu.val in mu.vals) {
      
      # Simulate raw individual-level data (raw = TRUE) with large N
      sim_dat_pop[[as.character(a.val)]][[type]][[as.character(mu.val)]] <-
        tryCatch(
          {
            sim.estd(
              a2.val  = a.val,
              d2.val  = d2.val,
              i2.val  = i2.val,
              c2.val  = c2.val,
              t2.val  = t2.val,
              mu.val  = mu.val,
              Nmz     = 5000,   # MZ twin pairs
              Ndz     = 5000,   # DZ twin pairs
              am.type = type,
              raw     = TRUE,   # Return individual-level data rather than summary stats
              exact   = FALSE 
            )
          },
          error = function(e) {
            message(sprintf(
              "Skipping: a2 = %s, type = %s, mu = %s | Error: %s",
              a.val, type, mu.val, e$message
            ))
            return(NULL)
          }
        )
    }
  }
}

# Model fitting: ETFD (iAM-ADE) ####
# ETFD uses 5 (twin 1, co-twin 2, sibling, partner 1, partner 2)
mean_MZ     <- rep("mean", 5)
mean_DZ     <- rep("mean", 5)
sel_vars    <- c("stw1", "tw1", "s1", "tw2", "stw2")    # ETFD phenotype columns

est_etfd_all <- NULL

for (a.val in a.vals) {
  for (type in am.type) {       
    for (mu.val in mu.vals) {
      
      # Retrieve simulated dataset for this parameter combination
      df <- sim_dat_pop[[as.character(a.val)]][[type]][[as.character(mu.val)]]$dat
      
      if (is.null(df)) next       # Skip if simulation failed for this cell
      
      # Split by zygosity (1 = MZ, 2 = DZ)
      mz <- df %>% filter(zyg == 1)
      dz <- df %>% filter(zyg == 2)
      
      # Choose direct or indirect assortative mating model based on what scenario was generated:
      # fit different NFTD models
      iAM_ADE_mod <- if (type == "dAM") {
        iam.date.fun(mz, dz, sel_vars, mean.MZ = mean_MZ, mean.DZ = mean_DZ, indirect.assortment = F, t1.free = F, t1s.free = F, name = "dAM_ADE")
      } else {
        iam.date.fun(mz, dz, sel_vars, mean.MZ = mean_MZ, mean.DZ = mean_DZ, indirect.assortment = T, t1.free = F, t1s.free = F, name = "iAM_ADE")
      }
       
      
      # Fit the CTD model and extract standardised variance components
      est_etfd_i <- tryCatch(
        {
          mxTryHard(iAM_ADE_mod, intervals = FALSE)$est_var$result[
            , c("varP_A", "varP_D", "varP_E")
          ] %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            rename(component = rowname, est = ".") %>%
            mutate(
              model = "iAM-ADE",
              a2         = a.val,   # Tag with generating parameter values
              assortment = type,    # so rows remain identifiable after rbind
              mu         = mu.val
            )
        },
        error = function(e) {
          message(sprintf(
            "Model fit failed: a2 = %s, type = %s, mu = %s | Error: %s",
            a.val, type, mu.val, e$message
          ))
          return(NULL)
        }
      )
      
      # Accumulate this iteration's estimates into the growing results frame
      est_etfd_all <- rbind(est_etfd_all, est_etfd_i)
      
    }
  }
}

# bias in genetic effects
est_am_bias <- est_etfd_all %>%
  # select standardized components
  filter(grepl("var", component),
         mu != 0) %>%
  mutate(
    d2 = .15, # we simulated dominance at .15
    component = case_match(
      component,
      "varP_E" ~ "E",
      "varP_T" ~ "T",
      "varP_D" ~ "D",
      "varP_C" ~ "C",
      "varP_E" ~ "E",
      "varP_D" ~ "D",
      "varP_C" ~ "C",
      "varP_N" ~ "N",
      "varP_N" ~ "N",
      "varP_A" ~ "A",
      "varP_A" ~ "A"
    )
  ) %>%
  dplyr::filter(component %in% c("A", "D")) %>%
  pivot_wider(names_from = c(, "component"), values_from = "est")


## Figure 2e ####
plot_d <- est_am_bias %>%
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  filter(mu != 0) %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  ggplot(aes(a2, A, color = mu)) +
  geom_point() +
  geom_vline(xintercept = .30, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed", color = "#B41F01") +
  scale_color_gradientn(colours = c(
    "#120246",
    "#050353",
    "#06175F",
    "#08306B",
    "#285C80",
    "#478495",
    "#68A7A9",
    "#89BDB5",
    "#AAD0C4",
    "#CCE3D7"
  )) +
  ylim(c(.2, 1)) +
  facet_grid(cols = vars(assortment), scales = "free") +
  theme_bw_classic(base_size = 10) +
  scale_x_continuous(breaks = seq(.30, 1, .15), limits = c(.30, 0.85)) +
  labs(
    x = expression(italic(h)^2),
    y = expression(hat(italic(h))[twin - ped]^2),
    color = expression(paste("Partner correlation"))
  ) +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )

## Figure 2f ####
plot_e <- est_am_bias %>%
  filter(assortment == "dAM" | assortment == "broad-sense hom." | assortment == "narrow-sense hom." | assortment == "mixed hom.") %>%
  mutate(assortment = factor(assortment, levels = c("dAM", "narrow-sense hom.", "broad-sense hom.", "mixed hom."))) %>%
  filter(mu != 0) %>%
  ggplot(aes(mu, D, color = a2)) +
  geom_point() +
  ylim(c(-0.01, 0.3)) +
  geom_hline(yintercept = .15, linetype = "dashed", color = "#FAA53C", linewidth = 1) +
  scale_color_viridis_c(option = "plasma", breaks = seq(.30, 1, .15)) +
  scale_x_continuous(limits = c(0, .95), breaks = seq(0, 1, .2)) +
  facet_grid(cols = vars(assortment), scales = "free") +
  theme_bw_classic(base_size = 10) +
  labs(
    x = expression(italic(mu)),
    y = expression(hat(italic(delta))[twin - ped]^2),
    color = expression(paste("Partner correlation"))
  ) +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )


## Figure 2 ####
Fig_2 <- 
  (plot_a / (plot_b + plot_c) / (plot_d + plot_e)) +
  plot_layout(
    heights = c(1.25, 1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "a"
  )

ggsave(
  filename = sprintf("%s/%s/00_FIG_2.pdf", wd_oa, wd_oa_figures),
  plot     = Fig_2,
  width    = 10.75,
  height   = 6,
  device   = cairo_pdf
)

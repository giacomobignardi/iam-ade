# Helper for iAM-ADE estimates plotting
# Written by G.B
# This functions are not likely not very helpfull (need customization)

## Shared lookup maps #### 
recode_type <- c(
  "r_MZ"          = "rmz",
  "r_DZ"          = "rdz",
  "r_FS"          = "rsib",
  "r_inlaw_MZ"    = "rmz-in-law",
  "r_inlaw_DZ"    = "rdz-in-law",
  "r_inlaw_FS"    = "rsib-in-law",
  "r_coinlaw_MZ"  = "rmz-co-in-law",
  "r_coinlaw_DZ"  = "rdz-co-in-law",
  "r_Partner_P"   = "ram",
  "r_Partner_S"   = "ram_s",  
  "r_inter_MZ"          = "rmz(adj)",
  "r_inter_DZ"          = "rdz(adj)",
  "r_inter_FS"          = "rsib(adj)",
  "r_inlaw_inter_MZ"    = "rmz-in-law(adj)",
  "r_inlaw_inter_DZ"    = "rdz-in-law(adj)",
  "r_inlaw_inter_FS"    = "rsib-in-law(adj)",
  "r_coinlaw_inter_MZ"  = "rmz-co-in-law(adj)",
  "r_coinlaw_inter_DZ"  = "rdz-co-in-law(adj)"

)

type_levels <- rev(c(
  "rmz", "rdz", "rsib",
  "rmz-in-law", "rdz-in-law", "rsib-in-law",
  "rmz-co-in-law", "rdz-co-in-law",
  "ram", "ram_s"
))

cor_levels <- c("iAM-ACE", "iAM-ADE", "dAM-ACE", "xdAM-ADE", "qdAM-ADE", "dAM-ADE", "Observed")

type_labels <- c(
  "ram_s"         = "spouse (sorting)",
  "ram"           = "spouse (focal)",
  "rmz"           = "MZ",
  "rdz"           = "DZ",
  "rsib"          = "sibling",
  "rmz-in-law"    = "MZ in-law",
  "rsib-in-law"   = "sibling|DZ in-law",
  "rmz-co-in-law" = "MZ co-in-law",
  "rdz-co-in-law" = "DZ co-in-law"
)

theme_bw_classic <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid.major   = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      strip.background   = element_blank(),
      legend.title.align = 0.5
    )
}

## Generic panel-building functions #### 
#' Build panel A (family correlations)
#'
#' @param fam_cor_long  Combined observed + implied correlations (long format)
#' @param implied_label Model label used for the ram_s annotation (e.g. "iAM-ACE")
#' @param implied_color Colour for the ram_s annotation text
#' @param fill_values   Named colour vector for scale_fill_manual
#' @param shape_values  Named shape vector for scale_shape_manual
#' @param exclude_types Character vector of type codes to drop (default: "rdz-in-law")

make_panel_a <- function(fam_cor_long,
                         implied_label  = "iAM-ADE",
                         implied_color  = "#E59F00",
                         fill_values    = c("iAM-ACE"  = "#E59F00",
                                            "dAM-ACE"  = "#009E73",
                                            "Observed" = "#0c72b2"),
                         shape_values   = c(24, 22, 21),
                         exclude_types  = "rdz-in-law",
                         trait_name     = "trait",
                         y_min          = -0.2,
                         y_max          =  1.0) {
  
  fam_cor_long %>%
    mutate(
      type        = factor(type, levels = type_levels),
      correlation = factor(correlation, levels = cor_levels)
    ) %>%
    filter(!type %in% exclude_types) %>%
    ggplot(aes(x = type, y = estimate, fill = correlation, shape = correlation)) +
    
    geom_vline(xintercept = 1.5, color = "gray", linetype = "dashed") +
    geom_hline(yintercept = 0,   color = "lightgray") +
    
    geom_linerange(
      aes(ymin = lbound, ymax = ubound),
      color = "gray", size = 0.75, position = position_dodge(0.5)
    ) +
    geom_point(size = 2.5, position = position_dodge(0.5)) +
    
    # Labels for all observed types except ram_s
    geom_text(
      data      = ~ filter(., correlation == "Observed", type != "ram_s"),
      aes(label = sprintf("%.2f\n[%.2f, %.2f]", estimate, lbound, ubound)),
      y = -Inf, vjust = -0.3, size = 3, lineheight = 0.8
    ) +
    
    # Label for ram_s from the implied model
    geom_text(
      data      = ~ filter(., correlation == implied_label, type == "ram_s"),
      aes(label = sprintf("%.2f\n[%.2f, %.2f]", estimate, lbound, ubound)),
      color = implied_color, y = -Inf, vjust = -0.3, size = 3, lineheight = 0.8
    ) +
    
    scale_y_continuous(limits = c(y_min, y_max), breaks = seq(y_min, y_max, 0.2)) +
    scale_shape_manual(name = "Correlation", values = shape_values) +
    scale_fill_manual( name = "Correlation", values = fill_values) +
    scale_x_discrete(labels = function(x) str_wrap(type_labels[x], width = 8)) +
    
    facet_wrap(~title, nrow = 1) +
    labs(y = paste0("Estimate (",trait_name,")"), x = "Relationship type") +
    theme_bw_classic() +
    theme()
}

#' Build panel B (variance components bar chart)
#'
#' @param est          Filtered & labelled variance-component data frame
#' @param param_levels Ordered factor levels for `parameter`
#' @param x_labels     Named vector of expressions for scale_x_discrete
#' @param fill_colors  Colour vector (length = number of parameters)

make_panel_b <- function(est,
                         param_levels,
                         x_labels,
                         fill_colors) {
  
  est %>%
    mutate(parameter = factor(parameter, levels = param_levels)) %>%
    ggplot(aes(x = parameter, y = estimate, fill = parameter)) +
    
    geom_bar(stat = "identity", position = position_dodge(),
             color = "black", width = 0.75) +
    
    geom_errorbar(aes(ymin = lbound, ymax = ubound),
                  position = position_dodge(), width = 0.3) +
    
    geom_text(aes(label = round(estimate, 2), y = ubound),
              position = position_dodge(width = 0.9),
              vjust = -0.5, size = 2.75) +
    
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete(labels = x_labels) +
    scale_fill_manual(values = fill_colors, guide = "none") +
    
    facet_wrap(~trait, nrow = 1, scales = "free") +
    labs(y = "Proportion of variance", x = "Parameter") +
    theme_bw_classic() +
    theme(legend.position = "none")
}

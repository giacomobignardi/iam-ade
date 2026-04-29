# Helper for iAM-ADE model
# Written by G.B
# Adapted from the iAM-ACE by Sunde et al. (2025, Nat Commun)
# load functions and models to fit and plot iAM-CATE and iAM-DATE

# Safe model fitting ####----------------------------------------------------
# Fit model with mxRun; fall back to mxTryHard if status code is non-zero.
# Returns a list with the fit object, convergence status, and a flag indicating
# whether mxTryHard was used.

run.model.safe <- function(model) {
  
  # Initial fit
  fit          <- mxRun(model, silent = TRUE)
  used_tryhard <- FALSE
  
  # If non-zero status, retry with mxTryHard
  if (fit$output$status$code != 0) {
    message("Model '", model@name, "' status ", fit$output$status$code, " — running mxTryHard")
    fit          <- mxTryHard(model, extraTries = 10, silent = TRUE)
    used_tryhard <- TRUE
  }
  
  # Return fit object and diagnostics
  list(
    fit          = fit,
    used_tryhard = used_tryhard,
    status_code  = fit$output$status$code,
    status_msg   = fit$output$status$statusMsg,
    model_name   = model@name
  )
}

# Helper: extract a fitted model by name ####--------------------------------
get_fit <- function(name) {
  idx <- iam_diagnostics %>% filter(model_name == name) %>% pull(index)
  iam_results[[idx]]$fit
}

# xiAM-DATE (sibling-shared preferences) ####--------------------------------

## Indirect assortment ####

# Full model (DATE: dominance + twin-shared environment)
xiam_date_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = T,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              # focal paths
                              d1.free = T, c1.free = F, t1.free = T,
                              # sorting paths
                              d1s.free = T, c1s.free = F, t1s.free = T,
                              # starting values
                              d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                              # sibling-shared preferences
                              cor.preferences = T, preferences.mz.dz = T,
                              name = "xiam_date"
)

# No twin-shared environment (ADE)
xiam_ade_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = T,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = T, c1.free = F, t1.free = F,
                             d1s.free = T, c1s.free = F, t1s.free = F,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             cor.preferences = T, preferences.mz.dz = T,
                             name = "xiam_ade"
)

# No dominance (ATE)
xiam_ate_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = T,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = F, t1.free = T,
                             d1s.free = F, c1s.free = F, t1s.free = T,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             cor.preferences = T, preferences.mz.dz = T,
                             name = "xiam_ate"
)

# Simplest model (AE)
xiam_ae_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = T,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = F, t1.free = F,
                            d1s.free = F, c1s.free = F, t1s.free = F,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            cor.preferences = T, preferences.mz.dz = T,
                            name = "xiam_ae"
)

## Direct assortment ####

# Full model (DATE)
xdam_date_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = F,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              d1.free = T, c1.free = F, t1.free = T,
                              d1s.free = T, c1s.free = F, t1s.free = T,
                              d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                              cor.preferences = T, preferences.mz.dz = T,
                              name = "xdam_date"
)

# No twin-shared environment (ADE)
xdam_ade_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = T, c1.free = F, t1.free = F,
                             d1s.free = T, c1s.free = F, t1s.free = F,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             cor.preferences = T, preferences.mz.dz = T,
                             name = "xdam_ade"
)

# No dominance (ATE)
xdam_ate_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = F, t1.free = T,
                             d1s.free = F, c1s.free = F, t1s.free = T,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             cor.preferences = T, preferences.mz.dz = T,
                             name = "xdam_ate"
)

# Simplest model (AE)
xdam_ae_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = F,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = F, t1.free = F,
                            d1s.free = F, c1s.free = F, t1s.free = F,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            cor.preferences = T, preferences.mz.dz = T,
                            name = "xdam_ae"
)

# Model list
xiam_date_list <- list(
  xiam_date_mod, xiam_ade_mod, xiam_ate_mod, xiam_ae_mod,
  xdam_date_mod, xdam_ade_mod, xdam_ate_mod, xdam_ae_mod
)


# qiAM-DATE (stratification) ####-------------------------------------------

## Indirect assortment ####

# Full model (DATE)
qiam_date_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = T,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              d1.free = T, c1.free = F, t1.free = T,
                              d1s.free = T, c1s.free = F, t1s.free = T,
                              d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                              qp.free = T,
                              name = "qiam_date"
)

# No twin-shared environment (ADE)
qiam_ade_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = T,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = T, c1.free = F, t1.free = F,
                             d1s.free = T, c1s.free = F, t1s.free = F,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             qp.free = T,
                             name = "qiam_ade"
)

# No dominance (ATE)
qiam_ate_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = T,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = F, t1.free = T,
                             d1s.free = F, c1s.free = F, t1s.free = T,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             qp.free = T,
                             name = "qiam_ate"
)

# Simplest model (AE)
qiam_ae_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = T,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = F, t1.free = F,
                            d1s.free = F, c1s.free = F, t1s.free = F,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            qp.free = T,
                            name = "qiam_ae"
)

## Direct assortment ####

# Full model (DATE)
qdam_date_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = F,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              d1.free = T, c1.free = F, t1.free = T,
                              d1s.free = T, c1s.free = F, t1s.free = T,
                              d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                              qp.free = T,
                              name = "qdam_date"
)

# No twin-shared environment (ADE)
qdam_ade_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = T, c1.free = F, t1.free = F,
                             d1s.free = T, c1s.free = F, t1s.free = F,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             qp.free = T,
                             name = "qdam_ade"
)

# No dominance (ATE)
qdam_ate_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = F, t1.free = T,
                             d1s.free = F, c1s.free = F, t1s.free = T,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             qp.free = T,
                             name = "qdam_ate"
)

# Simplest model (AE)
qdam_ae_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = F,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = F, t1.free = F,
                            d1s.free = F, c1s.free = F, t1s.free = F,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            qp.free = T,
                            name = "qdam_ae"
)

# Model list
qiam_date_list <- list(
  qiam_date_mod, qiam_ade_mod, qiam_ate_mod, qiam_ae_mod,
  qdam_date_mod, qdam_ade_mod, qdam_ate_mod, qdam_ae_mod
)


# iAM-DATE ####--------------------------------------------------------------

## Indirect assortment ####

# Full model (DATE)
iam_date_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = T,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = T, c1.free = F, t1.free = T,
                             d1s.free = T, c1s.free = F, t1s.free = T,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             name = "iam_date"
)

# No twin-shared environment (ADE)
iam_ade_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = T,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = T, c1.free = F, t1.free = F,
                            d1s.free = T, c1s.free = F, t1s.free = F,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            name = "iam_ade"
)

# No dominance (ATE)
iam_ate_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = T,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = F, t1.free = T,
                            d1s.free = F, c1s.free = F, t1s.free = T,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            name = "iam_ate"
)

# Simplest model (AE)
iam_ae_mod <- iam.date.fun(mz, dz, sel_vars,
                           indirect.assortment = T,
                           mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                           d1.free = F, c1.free = F, t1.free = F,
                           d1s.free = F, c1s.free = F, t1s.free = F,
                           d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                           name = "iam_ae"
)

## Direct assortment ####

# Full model (DATE)
dam_date_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = T, c1.free = F, t1.free = T,
                             d1s.free = T, c1s.free = F, t1s.free = T,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             name = "dam_date"
)

# No twin-shared environment (ADE)
dam_ade_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = F,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = T, c1.free = F, t1.free = F,
                            d1s.free = T, c1s.free = F, t1s.free = F,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            name = "dam_ade"
)

# No dominance (ATE)
dam_ate_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = F,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = F, t1.free = T,
                            d1s.free = F, c1s.free = F, t1s.free = T,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            name = "dam_ate"
)

# Simplest model (AE)
dam_ae_mod <- iam.date.fun(mz, dz, sel_vars,
                           indirect.assortment = F,
                           mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                           d1.free = F, c1.free = F, t1.free = F,
                           d1s.free = F, c1s.free = F, t1s.free = F,
                           d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                           name = "dam_ae"
)

# Model list
iam_date_list <- list(
  iam_date_mod, iam_ade_mod, iam_ate_mod, iam_ae_mod,
  dam_date_mod, dam_ade_mod, dam_ate_mod, dam_ae_mod
)


# DATE ####--------------------------------------------------------------
## Random assortment ####

# Full model (DATE)
date_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F, 
                             mu.free = F, mu.val = 0,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = T, c1.free = F, t1.free = T,
                             d1s.free = T, c1s.free = F, t1s.free = T,
                             d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                             name = "date"
)

# No twin-shared environment (ADE)
ade_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = F,
                            mu.free = F, mu.val = 0,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = T, c1.free = F, t1.free = F,
                            d1s.free = T, c1s.free = F, t1s.free = F,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            name = "ade"
)

# No dominance (ATE)
ate_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = F,
                            mu.free = F, mu.val = 0,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = F, t1.free = T,
                            d1s.free = F, c1s.free = F, t1s.free = T,
                            d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                            name = "ate"
)

# Simplest model (AE)
ae_mod <- iam.date.fun(mz, dz, sel_vars,
                           indirect.assortment = F,
                            mu.free = F, mu.val = 0,
                           mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                           d1.free = F, c1.free = F, t1.free = F,
                           d1s.free = F, c1s.free = F, t1s.free = F,
                           d1.val = .0, a1.val = .8, c1.val = .0, t1.val = .0, e1.val = .3,
                           name = "ae"
)

# Model list
date_list <- list(
  date_mod, ade_mod, ate_mod, ae_mod
)

# xiAM-CATE (sibling-shared preferences) ####--------------------------------

## Indirect assortment ####

# Full model (CATE: shared environment + twin-shared environment)
xiam_cate_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = T,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              d1.free = F, c1.free = T, t1.free = T,
                              d1s.free = F, c1s.free = T, t1s.free = T,
                              d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                              cor.preferences = T, preferences.mz.dz = T,
                              name = "xiam_cate"
)

# No twin-shared environment (ACE)
xiam_ace_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = T,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = T, t1.free = F,
                             d1s.free = F, c1s.free = T, t1s.free = F,
                             d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                             cor.preferences = T, preferences.mz.dz = T,
                             name = "xiam_ace"
)

# Note: ATE and AE models reuse xiam_ate_mod and xiam_ae_mod from section 2

## Direct assortment ####

# Full model (CATE)
xdam_cate_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = F,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              d1.free = F, c1.free = T, t1.free = T,
                              d1s.free = F, c1s.free = T, t1s.free = T,
                              d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                              cor.preferences = T, preferences.mz.dz = T,
                              name = "xdam_cate"
)

# No twin-shared environment (ACE)
xdam_ace_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = T, t1.free = F,
                             d1s.free = F, c1s.free = T, t1s.free = F,
                             d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                             cor.preferences = T, preferences.mz.dz = T,
                             name = "xdam_ace"
)

# Model list
xiam_cate_list <- list(
  xiam_cate_mod, xiam_ace_mod, xiam_ate_mod, xiam_ae_mod,
  xdam_cate_mod, xdam_ace_mod, xdam_ate_mod, xdam_ae_mod
)


# qiAM-CATE (stratification) ####-------------------------------------------

## Indirect assortment ####

# Full model (CATE)
qiam_cate_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = T,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              d1.free = F, c1.free = T, t1.free = T,
                              d1s.free = F, c1s.free = T, t1s.free = T,
                              d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                              qp.free = T,
                              name = "qiam_cate"
)

# No twin-shared environment (ACE)
qiam_ace_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = T,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = T, t1.free = F,
                             d1s.free = F, c1s.free = T, t1s.free = F,
                             d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                             qp.free = T,
                             name = "qiam_ace"
)

# Note: ATE and AE models reuse qiam_ate_mod and qiam_ae_mod from section 3

## Direct assortment ####

# Full model (CATE)
qdam_cate_mod <- iam.date.fun(mz, dz, sel_vars,
                              indirect.assortment = F,
                              mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                              d1.free = F, c1.free = T, t1.free = T,
                              d1s.free = F, c1s.free = T, t1s.free = T,
                              d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                              qp.free = T,
                              name = "qdam_cate"
)

# No twin-shared environment (ACE)
qdam_ace_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = T, t1.free = F,
                             d1s.free = F, c1s.free = T, t1s.free = F,
                             d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                             qp.free = T,
                             name = "qdam_ace"
)

# Model list
qiam_cate_list <- list(
  qiam_cate_mod, qiam_ace_mod, qiam_ate_mod, qiam_ae_mod,
  qdam_cate_mod, qdam_ace_mod, qdam_ate_mod, qdam_ae_mod
)


# iAM-CATE ####--------------------------------------------------------------

## Indirect assortment ####

# Full model (CATE)
iam_cate_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = T,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = T, t1.free = T,
                             d1s.free = F, c1s.free = T, t1s.free = T,
                             d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                             name = "iam_cate"
)

# No twin-shared environment (ACE)
iam_ace_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = T,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = T, t1.free = F,
                            d1s.free = F, c1s.free = T, t1s.free = F,
                            d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                            name = "iam_ace"
)

# Note: ATE and AE models reuse iam_ate_mod and iam_ae_mod from section 4

## Direct assortment ####

# Full model (CATE)
dam_cate_mod <- iam.date.fun(mz, dz, sel_vars,
                             indirect.assortment = F,
                             mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                             d1.free = F, c1.free = T, t1.free = T,
                             d1s.free = F, c1s.free = T, t1s.free = T,
                             d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                             name = "dam_cate"
)

# No twin-shared environment (ACE)
dam_ace_mod <- iam.date.fun(mz, dz, sel_vars,
                            indirect.assortment = F,
                            mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                            d1.free = F, c1.free = T, t1.free = F,
                            d1s.free = F, c1s.free = T, t1s.free = F,
                            d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                            name = "dam_ace"
)

# Model list
iam_cate_list <- list(
  iam_cate_mod, iam_ace_mod, iam_ate_mod, iam_ae_mod,
  dam_cate_mod, dam_ace_mod, dam_ate_mod, dam_ae_mod
)

# CATE ####--------------------------------------------------------------
## Radnom assortment ####

# Full model (CATE)
cate_mod <- iam.date.fun(mz, dz, sel_vars,
                         indirect.assortment = F,
                         mu.free = F, mu.val = 0,
                         mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                         d1.free = F, c1.free = T, t1.free = T,
                         d1s.free = F, c1s.free = T, t1s.free = T,
                         d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                         name = "cate"
)

# No twin-shared environment (ACE)
ace_mod <- iam.date.fun(mz, dz, sel_vars,
                        indirect.assortment = F,
                        mu.free = F, mu.val = 0,
                        mean.MZ = mean_MZ, mean.DZ = mean_DZ,
                        d1.free = F, c1.free = T, t1.free = F,
                        d1s.free = F, c1s.free = T, t1s.free = F,
                        d1.val = .0, a1.val = .6, c1.val = .2, t1.val = .0, e1.val = .3,
                        name = "ace"
)

# Model list
cate_list <- list(
  cate_mod, ace_mod, ate_mod, ae_mod
)
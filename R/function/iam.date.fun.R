# iAM-DATE function
# Written by Hans Fredrik Sunde (H.F.S) & Giacomo Bignardi (G.B)
# Adapted from the iAM-ACE by Sunde et al. (2025, Nat Commun)
# Partner of twins model (MZ and DZ twins, their spouses and one sibling)
# With both iAM and dAM (indirect and direct assortment)
# Version 04.2026

message("Version 04.2026")
mxOption(key = "Default optimizer", value = "NPSOL")

iam.date.fun <- function(
    
  # Data
  mzData  = NULL,
  dzData  = NULL,
  selVars = c("S1_p", "T1_p", "FS1_p", "T2_p", "S2_p"), # exact order, names can be changed
  name    = "iAM_DATE",
  
  # Focal phenotype
  a1.free = TRUE,  # Freely estimate a?
  c1.free = FALSE, # Freely estimate c?
  t1.free = TRUE,  # Freely estimate t?
  d1.free = TRUE,  # Freely estimate d?
  n1.free = FALSE, # Freely estimate n?
  e1.free = TRUE,  # Freely estimate e?
  a1.val  = .7,    # Starting value for a
  c1.val  = .0,    # Starting value for c
  t1.val  = .0,    # Starting value for t
  d1.val  = .0,    # Starting value for d
  n1.val  = .0,    # Starting value for n
  e1.val  = .7,    # Starting value for e
  

  # Assortative mating
  mu.free             = TRUE,  # Freely estimate assortment strength (mu)?
  mu.val              = .3,    # Starting value for mu
  U.f                 = 1,     # Intergenerational equilibrium (0 = first generation, 1 = equilibrium)
  indirect.assortment = FALSE, # If FALSE, sorting factor is equated to focal phenotype
  VarS_constraint     = TRUE,  # Constrain Var(S) == Var(P); alternative: fix one path (e.g. e1s)
  
  # Sorting factor
  # If fixed and NULL: tied to corresponding focal phenotype path
  # If fixed and non-NULL: constrained to that value
  a1s.free = TRUE,  # Freely estimate scaling factor for VA? (genetic homogamy, additive)
  c1s.free = FALSE, # Freely estimate scaling factor for VC? (social homogamy, sibling-shared)
  t1s.free = TRUE,  # Freely estimate scaling factor for VT? (social homogamy, twin-shared)
  d1s.free = TRUE,  # Freely estimate scaling factor for VD? (genetic homogamy, dominance)
  n1s.free = FALSE, # Freely estimate scaling factor for VN? (genetic homogamy, epistasis)
  e1s.free = TRUE,  # Freely estimate scaling factor for VE? (idiosyncratic homogamy; must be fixed unless others are fixed too)
  a1s.val  = NULL,
  c1s.val  = NULL,
  t1s.val  = NULL,
  d1s.val  = NULL,
  n1s.val  = NULL,
  e1s.val  = NULL,
  
  # Stratification
  qp.free           = FALSE, # Freely estimate VQ? (Stratification)
  
  # Mate preferences
  cor.preferences   = FALSE, # Freely estimate sibling-correlated mate preferences?
  preferences.mz.dz = FALSE, # Allow preference-correlation to differ across zygosity?
  
  # Miscellaneous
  parental_rGE = 0,              # Gene-environment correlation in parental generation (not yet implemented)
  ni           = 2,              # AxA interaction loci (rAxA = 1 / 2^ni; default = 2-locus AxA) # USE WITH CAUTION: assumes no increase in epistatic variance under assortative mating
  mean.MZ      = c("mean", "mean", "mean", "mean", "mean"), # Mean labels for MZ group
  mean.DZ      = c("mean", "mean", "mean", "mean", "mean"), # Mean labels for DZ group
  fitFunction  = "ML",           # Fit function (currently only ML is supported)
  
  # Measurement error
  # USE WITH CAUTION
  # Focal phenotype error
  r1.free = FALSE,           # Estimate error?
  ryy     = 1,               # test-retest reliability
  r1.val  = sqrt(1 - ryy)    # Standardised error/intra-individual path
  
) {
  
  timestamp()
  require(OpenMx)
  
  ntv <- length(selVars)
  
  # Input checks ---------------------------------------------------------------
  
  if (!cor.preferences) preferences.mz.dz <- FALSE # preferences.mz.dz can only be free if cor.preferences is included
  if (cor.preferences & qp.free) warning("Social stratification and correlated preferences cannot be estimated simultaneously. Expect errors and untrustworthy results.")
  if (a1s.free & c1s.free & t1s.free & e1s.free & !VarS_constraint) stop(
    "All paths to the sorting factor are free but VarS_constraint = FALSE.\n",
    "Fix at least one path (e.g. e1s.free = FALSE) or set VarS_constraint = TRUE."
  )
  
  # Resolve sorting factor labels and starting values --------------------------
  # If a path is fixed with no user-supplied value, it is tied to the corresponding focal phenotype path
  
  if (indirect.assortment) {
    if (!a1s.free & is.null(a1s.val)) { lab_a1s <- "est_a1";  a1s.free <- a1.free; a1s.val <- a1.val } else { lab_a1s <- "est_a1s"; if (is.null(a1s.val)) a1s.val <- a1.val }
    if (!c1s.free & is.null(c1s.val)) { lab_c1s <- "est_c1";  c1s.free <- c1.free; c1s.val <- c1.val } else { lab_c1s <- "est_c1s"; if (is.null(c1s.val)) c1s.val <- c1.val }
    if (!t1s.free & is.null(t1s.val)) { lab_t1s <- "est_t1";  t1s.free <- t1.free; t1s.val <- t1.val } else { lab_t1s <- "est_t1s"; if (is.null(t1s.val)) t1s.val <- t1.val }
    if (!d1s.free & is.null(d1s.val)) { lab_d1s <- "est_d1";  d1s.free <- d1.free; d1s.val <- d1.val } else { lab_d1s <- "est_d1s"; if (is.null(d1s.val)) d1s.val <- d1.val }
    if (!n1s.free & is.null(n1s.val)) { lab_n1s <- "est_n1";  n1s.free <- n1.free; n1s.val <- n1.val } else { lab_n1s <- "est_n1s"; if (is.null(n1s.val)) n1s.val <- n1.val }
    if (!e1s.free & is.null(e1s.val)) { lab_e1s <- "est_e1";  e1s.free <- e1.free; e1s.val <- e1.val } else { lab_e1s <- "est_e1s"; if (is.null(e1s.val)) e1s.val <- e1.val }
  } else {
    # Direct assortment: all sorting factor paths tied to focal phenotype
    lab_a1s <- "est_a1"; a1s.free <- a1.free; a1s.val <- a1.val
    lab_c1s <- "est_c1"; c1s.free <- c1.free; c1s.val <- c1.val
    lab_t1s <- "est_t1"; t1s.free <- t1.free; t1s.val <- t1.val
    lab_d1s <- "est_d1"; d1s.free <- d1.free; d1s.val <- d1.val
    lab_n1s <- "est_n1"; n1s.free <- n1.free; n1s.val <- n1.val
    lab_e1s <- "est_e1"; e1s.free <- e1.free; e1s.val <- e1.val
  }
  
  # Free parameters ------------------------------------------------------------
  
  paths <- list(
    
    # Fixed scalars
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = U.f, labels = "U1",  name = "U"),  # Intergenerational equilibrium
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = ni,  labels = "ni",  name = "nI"), # AxA interaction level
    
    # Focal phenotype (ACDTE)
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = a1.free, values = a1.val, labels = "est_a1", name = "a1"), # Additive genetic
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = c1.free, values = c1.val, labels = "est_c1", name = "c1"), # Sibling shared environment
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = t1.free, values = t1.val, labels = "est_t1", name = "t1"), # Twin shared environment
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = d1.free, values = d1.val, labels = "est_d1", name = "d1"), # Dominance genetic
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = n1.free, values = n1.val, labels = "est_n1", name = "n1"), # Epistatic genetic
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = e1.free, values = e1.val, labels = "est_e1", name = "e1"), # Unique environment
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE,   values = r1.val, labels = "est_r1", name = "r1"), # Intra-individual variance
    
    # Sorting factor (ACDTE scaling)
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = a1s.free, values = a1s.val, labels = lab_a1s, name = "a1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = c1s.free, values = c1s.val, labels = lab_c1s, name = "c1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = t1s.free, values = t1s.val, labels = lab_t1s, name = "t1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = d1s.free, values = d1s.val, labels = lab_d1s, name = "d1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = n1s.free, values = n1s.val, labels = lab_n1s, name = "n1s"),
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = e1s.free, values = e1s.val, labels = lab_e1s, name = "e1s"),
    
    # Sorting factor variance and covariance with focal phenotype
    mxAlgebra(expression = a1s^2 + c1s^2 + e1s^2 + 2*a1s*w*c1s + t1s^2 + d1s^2 + n1s^2,                           name = "VarS"),
    mxAlgebra(expression = a1*a1s + c1*c1s + e1*e1s + c1*w*a1s + a1*w*c1s + t1*t1s + d1*d1s + n1*n1s,             name = "covPS"),
    mxAlgebra(expression = mu * (a1s + c1s*w)^2,                                                                  name = "partner_rG"),
    mxAlgebra(expression = mu * (n1s)^2,                                                                          name = "partner_rN"),
    
    # Assortment strength (mu_param is the raw parameter; mu is rescaled to copath)
    mxMatrix(type = "Lower", nrow = 1, ncol = 1, free = mu.free, values = mu.val, labels = "mu1", name = "mu_param"),
    mxAlgebra(expression = mu_param / VarS^2, name = "mu"),
    
    # Social stratification
    mxMatrix(type = "Lower", nrow = 1, ncol = 1, free = qp.free,           labels = "q1", name = "q"),
    
    # Sibling-correlated mate preferences
    mxMatrix(type = "Lower", nrow = 1, ncol = 1, free = cor.preferences,   values = 0, labels = "z1", name = "Z_dz"), # Sibling-correlated preferences
    mxMatrix(type = "Lower", nrow = 1, ncol = 1, free = preferences.mz.dz, values = 0, labels = "z2", name = "Z_mz"), # MZ-specific preferences
    mxAlgebra(expression = Z_dz^2,       name = "Zdz"),
    mxAlgebra(expression = Zdz + Z_mz^2, name = "Zmz"),
    
    # Gene-environment covariance (parental generation)
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = parental_rGE, labels = "w_est",  name = "w"),
    
    # Fixed measurement error/intra-individual path (scaled)
    mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = r1.val,      labels = "val_zr1", name = "zr1")
  )
  
  # Constraints ----------------------------------------------------------------
  
  # Scale error path to trait variance (only applied if measurement error is non-zero)
  r1_cons   <- if (r1.free == T)    mxConstraint(expression = r1   == zry * sqrt(var1), name = "r1")   else NULL
  # Constrain Var(S) == Var(P); note var_bw1 == var1 when ryy = 1 (default)
  VarS_cons <- if (VarS_constraint) mxConstraint(expression = VarS == var1,             name = "VarS_constraint") else NULL
  
  # Variance components --------------------------------------------------------
  
  variances <- list(
    
    # Focal phenotype
    mxAlgebra(expression = a1^2,          name = "VA1"),
    mxAlgebra(expression = c1^2,          name = "VC1"),
    mxAlgebra(expression = t1^2,          name = "VT1"),
    mxAlgebra(expression = d1^2,          name = "VD1"),
    mxAlgebra(expression = n1^2,          name = "VN1"),
    mxAlgebra(expression = 2*c1*w*a1,     name = "VrGE1"),
    mxAlgebra(expression = e1^2,          name = "VE1"),
    mxAlgebra(expression = q^2,           name = "VQ1"),
    mxAlgebra(expression = r1^2,          name = "VR1"),
    
    # Sorting factor
    mxAlgebra(expression = a1s^2,         name = "VA1_s"),
    mxAlgebra(expression = c1s^2,         name = "VC1_s"),
    mxAlgebra(expression = t1s^2,         name = "VT1_s"),
    mxAlgebra(expression = d1s^2,         name = "VD1_s"),
    mxAlgebra(expression = n1s^2,         name = "VN1_s"),
    mxAlgebra(expression = 2*c1s*w*a1s,   name = "VrGE1_s"),
    mxAlgebra(expression = e1s^2,         name = "VE1_s")
  )
  
  # Expected variances ---------------------------------------------------------
  
  var1    <- mxAlgebra(expression = VA1 + VC1 + VrGE1 + VE1 + VT1 + VD1 + VN1 + VQ1 + VR1, name = "var1")    # Total variance
  var_bw1 <- mxAlgebra(expression = VA1 + VC1 + VrGE1 + VE1 + VT1 + VD1 + VN1 + VQ1,       name = "var_bw1") # Variance excluding intra-individual
  
  # Sibling genetic correlations -----------------------------------------------
  
  FS_rA <- mxAlgebra(expression = (1 + U*partner_rG) / 2,      name = "FS_rA") # Additive genetic correlation between siblings
  FS_rD <- mxAlgebra(expression = 1/4,                         name = "FS_rD") # Dominance genetic correlation between siblings
  FS_rN <- mxAlgebra(expression = (1 + U*partner_rN) / (2^nI), name = "FS_rN") # Epistatic genetic correlation between siblings
  
  # Expected covariances -------------------------------------------------------
  
  # Twin/sibling pairs
  cov_MZ <- mxAlgebra(expression =       VA1 + VC1 + VrGE1 +       VD1 +       VN1 + VT1 + VQ1, name = "cov_MZ") # MZ twins
  cov_DZ <- mxAlgebra(expression = FS_rA*VA1 + VC1 + VrGE1 + FS_rD*VD1 + FS_rN*VN1 + VT1 + VQ1, name = "cov_DZ") # DZ twins
  cov_FS <- mxAlgebra(expression = FS_rA*VA1 + VC1 + VrGE1 + FS_rD*VD1 + FS_rN*VN1       + VQ1, name = "cov_FS") # Full siblings
  
  # In-law pairs
  cov_In_MZ <- mxAlgebra(expression = covPS*mu * (      a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s +       d1*d1s +       n1*n1s + t1*t1s) + VQ1, name = "cov_In_MZ") # MZ in-law
  cov_In_DZ <- mxAlgebra(expression = covPS*mu * (FS_rA*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + FS_rD*d1*d1s + FS_rN*n1*n1s + t1*t1s) + VQ1, name = "cov_In_DZ") # DZ in-law
  cov_In_FS <- mxAlgebra(expression = covPS*mu * (FS_rA*a1*a1s + c1*c1s + c1s*w*a1 + c1*w*a1s + FS_rD*d1*d1s + FS_rN*n1*n1s         ) + VQ1, name = "cov_In_FS") # Sibling in-law
  
  # Co-in-law pairs
  cov_CIn_MZ <- mxAlgebra(expression = covPS^2 * (mu^2 * (      a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2 +       d1s^2 +       n1s^2) + Zmz) + VQ1, name = "cov_CIn_MZ") # MZ co-in-law
  cov_CIn_DZ <- mxAlgebra(expression = covPS^2 * (mu^2 * (FS_rA*a1s^2 + c1s^2 + 2*c1s*w*a1s + t1s^2 + FS_rD*d1s^2 + FS_rN*n1s^2) + Zdz) + VQ1, name = "cov_CIn_DZ") # DZ co-in-law
  
  # Partners
  cov_Mate <- mxAlgebra(expression = covPS*mu*covPS + VQ1, name = "cov_Mate")
  
  # Expected covariance matrices -----------------------------------------------
  # Variable order: S1 (spouse 1), T1 (twin 1), FS (sibling), T2 (twin 2), S2 (spouse 2)
  
  expCov_MZ <- mxAlgebra(
    expression = rbind(
      cbind(var1,       cov_Mate,  cov_In_FS, cov_In_MZ, cov_CIn_MZ),
      cbind(cov_Mate,   var1,      cov_FS,    cov_MZ,    cov_In_MZ),
      cbind(cov_In_FS,  cov_FS,    var1,      cov_FS,    cov_In_FS),
      cbind(cov_In_MZ,  cov_MZ,    cov_FS,    var1,      cov_Mate),
      cbind(cov_CIn_MZ, cov_In_MZ, cov_In_FS, cov_Mate,  var1)
    ),
    name = "expCov_MZ"
  )
  
  expCov_DZ <- mxAlgebra(
    expression = rbind(
      cbind(var1,       cov_Mate,  cov_In_FS, cov_In_DZ, cov_CIn_DZ),
      cbind(cov_Mate,   var1,      cov_FS,    cov_DZ,    cov_In_DZ),
      cbind(cov_In_FS,  cov_FS,    var1,      cov_FS,    cov_In_FS),
      cbind(cov_In_DZ,  cov_DZ,    cov_FS,    var1,      cov_Mate),
      cbind(cov_CIn_DZ, cov_In_DZ, cov_In_FS, cov_Mate,  var1)
    ),
    name = "expCov_DZ"
  )
  
  # Means and expectations -----------------------------------------------------
  
  mean_MZ <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = 0, labels = mean.MZ, name = "mean_MZ")
  mean_DZ <- mxMatrix(type = "Full", nrow = 1, ncol = 5, free = TRUE, values = 0, labels = mean.DZ, name = "mean_DZ")
  exp_MZ  <- mxExpectationNormal(covariance = "expCov_MZ", means = "mean_MZ", dimnames = selVars)
  exp_DZ  <- mxExpectationNormal(covariance = "expCov_DZ", means = "mean_DZ", dimnames = selVars)
  
  # Data -----------------------------------------------------------------------
  
  data_MZ <- mxData(observed = mzData, type = "raw")
  data_DZ <- mxData(observed = dzData, type = "raw")
  
  # Object lists ---------------------------------------------------------------
  
  common_objects <- list(var1, var_bw1, cov_Mate, FS_rA, FS_rD, FS_rN)
  exp_cov_list   <- list(expCov_MZ, expCov_DZ)
  algebras_MZ    <- list(mean_MZ, cov_MZ, cov_FS, cov_In_FS, cov_In_MZ, cov_CIn_MZ)
  algebras_DZ    <- list(mean_DZ, cov_DZ, cov_FS, cov_In_FS, cov_In_DZ, cov_CIn_DZ)
  algebras_all   <- list(algebras_MZ, algebras_DZ)
  
  # Group models ---------------------------------------------------------------
  
  fitFun   <- mxFitFunctionML()
  multi    <- mxFitFunctionMultigroup(c("MZ", "DZ"))
  model_MZ <- mxModel(paths, variances, common_objects, algebras_MZ, expCov_MZ, exp_MZ, data_MZ, fitFun, name = "MZ")
  model_DZ <- mxModel(paths, variances, common_objects, algebras_DZ, expCov_DZ, exp_DZ, data_DZ, fitFun, name = "DZ")
  model_list <- list(model_MZ, model_DZ)
  
  # Results algebras and confidence intervals ----------------------------------
  
  # Variance components (proportions of total and between-individual variance)
  est_var_Col <<- c(
    "varP", "varP_bw", "varS",
    "varP_A",    "varP_C",    "varP_T",    "varP_D",    "varP_N",    "varP_E",    "varP_Q",    "varP_Ew",
    "varP_bw_A", "varP_bw_C", "varP_bw_T", "varP_bw_D", "varP_bw_N", "varP_bw_E", "varP_bw_Q",
    "varS_A",    "varS_C",    "varS_T",    "varS_D",    "varS_N",    "varS_E"
  )
  
  est_var <- mxAlgebra(
    expression = cbind(
      var1, var_bw1, VarS,
      VA1/var1,    VC1/var1,    VT1/var1,    VD1/var1,    VN1/var1,    VE1/var1,    VQ1/var1,    VR1/var1,
      VA1/var_bw1, VC1/var_bw1, VT1/var_bw1, VD1/var_bw1, VN1/var_bw1, VE1/var_bw1, VQ1/var_bw1,
      VA1_s/VarS,  VC1_s/VarS,  VT1_s/VarS,  VD1_s/VarS,  VN1_s/VarS,  VE1_s/VarS
    ),
    name = "est_var", dimnames = list("est", est_var_Col)
  )
  
  # Correlations
  est_cov_Col <<- c(
    "r_MZ",             "r_DZ",             "r_FS",
    "r_bw_MZ",          "r_bw_DZ",          "r_bw_FS",
    "r_inlaw_MZ",       "r_inlaw_DZ",       "r_inlaw_FS",
    "r_inlaw_bw_MZ",    "r_inlaw_bw_DZ",    "r_inlaw_bw_FS",
    "r_coinlaw_MZ",     "r_coinlaw_DZ",
    "r_CIn_bw_MZ",      "r_CIn_bw_DZ",
    "r_Partner_P",      "r_Partner_S", 
    "rA_Partner",       "rN_Partner", 
    "rA_FS",            "rN_FS",
    "X_mz",   "X_dz"
  )
  
  est_cov <- mxAlgebra(
    expression = cbind(
      cov_MZ/var1,        cov_DZ/var1,       cov_FS/var1,
      cov_MZ/var_bw1,     cov_DZ/var_bw1,    cov_FS/var_bw1,
      cov_In_MZ/var1,     cov_In_DZ/var1,    cov_In_FS/var1,
      cov_In_MZ/var_bw1,  cov_In_DZ/var_bw1, cov_In_FS/var_bw1,
      cov_CIn_MZ/var1,    cov_CIn_DZ/var1,
      cov_CIn_MZ/var_bw1, cov_CIn_DZ/var_bw1,
      cov_Mate/var1,  mu*VarS, 
      partner_rG,     partner_rN, 
      FS_rA,          FS_rN,
      Zmz,Zdz
    ),
    name = "est_cov", dimnames = list("est", est_cov_Col)
  )
  
  # Path coefficients
  est_path_Col <<- c(
    "mu",
    "a1s", "c1s", "t1s", "d1s", "n1s", "e1s",
    "a1",  "c1",  "t1",  "d1",  "n1",  "e1",  "q",  "ew1"
  )
  
  est_path <- mxAlgebra(
    expression = cbind(
      mu_param,
      a1s, c1s, t1s, d1s, n1s, e1s,
      a1,  c1,  t1,  d1,  n1,  e1,  q1,  r1
    ),
    name = "est_path", dimnames = list("est", est_path_Col)
  )
  
  CI       <- mxCI(c("est_path", "est_var", "est_cov"))
  CI_names <<- c(est_path_Col, est_var_Col, est_cov_Col)
  
  results_objects <- list(est_path, est_var, est_cov)
  
  # Assemble model -------------------------------------------------------------
  
  model_iamDATE <- mxModel(
    paths, multi, VarS_cons, common_objects, variances,
    algebras_all, model_list, exp_cov_list, CI, results_objects,
    name = name
  )
  
}
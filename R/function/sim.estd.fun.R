# Simulate extended families of twin data
# Written by G.B
# Function to generate twin, sibling, and spouses correlations 
# under direct and indirect assortative mating, preferential mating, stratification, and measurement error

sim.estd <-  function(
    a2.val = .40,    # additive genetics
    d2.val = .00,    # dominance genetics
    i2.val = .10,    # other non-additive genetic effects beyond dominance - here interaction deviation
    c2.val = .05,    # shared environement
    t2.val = .00,    # twin-shared environement
    mu.val = .00,    # assortment between partners co-path coefficient
    xi = 0,          # cross-trait assortative mating on sibling-shared preferences (product of mu_ss and r_aes)
    q2.val = 0,      # stratification
    ryy = 1,         # test-retest reliability (between individual variance)
    am.type = "dAM", # type of assortment
    Nmz = 6,         # minimum sample to generate output with exact set to TRUE
    Ndz = 6,         # minimum sample to generate output with exact set to TRUE
    exact = T,       # if TRUE the parameter recovery should be exact (sample parameter, not estimate from population sample)
    raw   = F        # if FALSE the raw generated data are not provided
){
  ## starting values ####
  # as raw ML implies division by N not (N-1)
  cor_Nmz <- (Nmz - 1) / Nmz
  cor_Ndz <- (Ndz - 1) / Ndz
  N <- Nmz + Ndz
  Nmem <- 5 # 2 partners + 2 twins + 1 full-sib
  Nmem1 <- Nmem + 1
  exact <- exact # if TRUE the parameter recovery should be exact (sample parameter, not estimate from population sample)
  
  # set squared paths coefficients 
  a2 <- a2.val*(ryy)
  c2 <- c2.val*(ryy)
  t2 <- t2.val*(ryy)
  d2 <- d2.val*(ryy)
  i2 <- i2.val*(ryy)
  q2 <- q2.val*(ryy)
  ew2 <- 1-ryy
  
  # focal phenotype
  a1 <- sqrt(a2)     # additive genetics
  c1 <- sqrt(c2)     # shared environment (always assumed to be equal to 0)
  t1 <- sqrt(t2)     # twin environment
  d1 <- sqrt(d2)     # dominance genetics
  i1 <- sqrt(i2)     # interaction deviation
  q1 <- sqrt(q2)     # stratification
  ew1 <- sqrt(ew2) # unique environment (intra)
  eb1 <- sqrt(1 - a2 - d2 - i2 - t2 - ew2) # unique environment (inter)
  
  # build variance of the focal trait (relevant to the sorting trati)
  Vp_bt <- a1^2 + c1^2 + t1^2 + d1^2 + i1^2 + eb1^2 + q1^2            # phenotypic between individual variance focal phenotype (trait) 
  
  if (am.type == "dAM") {
    # sorting phenotype
    a1s  <- sqrt(a2)  # narrow-sense hom.
    c1s  <- sqrt(c2)  # social hom. (shared)
    t1s  <- sqrt(t2)  # social hom. (twin-shared)
    d1s  <- sqrt(d2)
    i1s  <- sqrt(i2)
    eb1s <- eb1       # idyosincratic hom.
  } else if (am.type == "narrow-sense hom.") {
    # sorting phenotype (narrow-sense hom.)
    a1s  <- sqrt(1)*sqrt(Vp_bt) # purely correlate via A
    c1s  <- log(1.0)
    t1s  <- log(1.0)
    d1s  <- log(1.0)
    i1s  <- log(1.0)
    eb1s <- log(1.0)
  } else if (am.type == "broad-sense hom.") {
    # sorting phenotype (broad-sense hom.)
    a1s  <- sqrt(2 / 3)*sqrt(Vp_bt) # correlate via A
    c1s  <- log(1.0)
    t1s  <- log(1.0)
    d1s  <- sqrt(1 / 3)*sqrt(Vp_bt) # correlate via D
    i1s  <- log(1.0)
    eb1s <- log(1.0)
  } else if (am.type == "twin-social hom.") {
    # sorting phenotype (twin-shared hom.)
    a1s  <- log(1.0)
    c1s  <- log(1.0)
    t1s  <- sqrt(1)*sqrt(Vp_bt) # purely correlate via T
    d1s  <- log(1.0)
    i1s  <- log(1.0)
    eb1s <- log(1.0)
  } else if (am.type == "social hom.") {
    # sorting phenotype (twin-shared hom.)
    a1s  <- log(1.0)
    c1s  <- sqrt(1)*sqrt(Vp_bt) # purely correlate via C
    t1s  <- log(1.0)
    d1s  <- log(1.0)
    i1s  <- log(1.0)
    eb1s <- log(1.0)
  } else if (am.type == "mixed hom.") {
    # sorting phenotype (mixed hom.)
    a1s  <- sqrt(4 / 9)*sqrt(Vp_bt) # correlate via A
    c1s  <- log(1.0)
    t1s  <- log(1.0) 
    d1s  <- sqrt(2 / 9)*sqrt(Vp_bt) # correlate via D
    i1s  <- log(1.0)
    eb1s <- sqrt(1 / 3)*sqrt(Vp_bt) # correlate via E
  } else if (am.type == "idio. hom.") {
    # sorting phenotype (idio. hom. homogamy)
    a1s  <- log(1.0)
    c1s  <- log(1.0)
    t1s  <- log(1.0)
    d1s  <- log(1.0)
    i1s  <- log(1.0)
    eb1s <- sqrt(1)*sqrt(Vp_bt) # purely correlate via E
  }
  
  ## generative model ####
  # build variances and co-pahts
  Vp <- a1^2 + c1^2 + t1^2 + d1^2 + i1^2 + eb1^2 + ew1^2 + q1^2       # phenotypic variance focal phenotype (trait)
  Vs <- a1s^2 + c1s^2 + t1s^2 + d1s^2 + i1s^2 + eb1s^2                # phenotypic variance sorting phenotype (sorting factor)
  mu <- (mu.val*(ryy^2)) / (Vs^2)                                     # mu value given is standardised - needs to be scaled if ryy != 1
  
  # randomly generate values for r_aes and mu_ss given their product, xi
  r_aes <- runif(1, min = xi, max = 1)
  mu_ss2 <- xi / r_aes
  
  # partner genetic correlations
  partner_rG <- mu * a1s^2
  
  # DZ and sib correlations
  fs_rA <- (1 + partner_rG) / 2 # A correlation
  fs_rD <- .25 # D correlation
  fs_rN <- .0 # N high-order epistasis boundary
  
  # build the MZ and DZ covariance matrices based on the parameter values given above
  Cps <- a1 * a1s + c1 * c1s + t1 * t1s + d1 * d1s + i1 * i1s + eb1 * eb1s # covariance between focal phenotype and sorting factor
  Cpar <- (Cps * mu * Cps) + q1^2 # covariance between partners
  Cmz <- 1 * a1^2 + c1^2 + 1 * d1^2 + 1 * i1^2 + t1^2 + q1^2 # mz covariance
  Cdz <- fs_rA * a1^2 + c1^2 + fs_rD * d1^2 + fs_rN * i1^2 + t1^2 + q1^2 # dz covariance
  Cfs <- fs_rA * a1^2 + c1^2 + fs_rD * d1^2 + fs_rN * i1^2 + q1^2 # fs covariance
  Cmzil <- mu * Cps * (1 * a1 * a1s + c1 * c1s + 1 * d1 * d1s + 1 * i1 * i1s + t1 * t1s) + q1^2 # mz in-law
  Cdzil <- mu * Cps * (fs_rA * a1 * a1s + c1 * c1s + fs_rD * d1 * d1s + fs_rN * i1 * i1s + t1 * t1s) + q1^2 # dz in-law
  Cfsil <- mu * Cps * (fs_rA * a1 * a1s + c1 * c1s + fs_rD * d1 * d1s + fs_rN * i1 * i1s) + q1^2 # fs in-law
  Cmzcil <- Cps^2 * (mu^2 * (1 * a1s^2 + c1s^2 + 1 * d1s^2 + 1 * i1s^2 + t1s^2) + mu_ss2*r_aes) + q1^2 # mz co-in-law
  Cdzcil <- Cps^2 * (mu^2 * (fs_rA * a1s^2 + c1s^2 + fs_rD * d1s^2 + fs_rN * i1s^2 + t1s^2) +mu_ss2*r_aes) + q1^2 # dz co-in-law
  
  # variance components
  # focal phenotype
  SA1  <- a1^2 / Vp
  SC1  <- c1^2 / Vp
  ST1  <- t1^2 / Vp
  SD1  <- d1^2 / Vp
  SI1  <- i1^2 / Vp
  SEW1 <- ew1^2 / Vp
  SEB1 <- eb1^2 / Vp
  SE1  <- (ew1^2 + eb1^2) / Vp
  
  # MZ covariance matrix
  # p1, tw1, s1, tw2, p2
  Smz <- matrix(
    c(
      Vp,     Cpar,  Cfsil, Cmzil, Cmzcil,
      Cpar,   Vp,    Cfs,   Cmz,   Cmzil,
      Cfsil,  Cfs,   Vp,    Cfs,   Cfsil,
      Cmzil,  Cmz,   Cfs,   Vp,    Cpar,
      Cmzcil, Cmzil, Cfsil, Cpar,  Vp
    ),
    5, 5,
    byrow = T
  )
  
  # DZ covariance matrix
  Sdz <- matrix(
    c(
      Vp,     Cpar,  Cfsil, Cdzil, Cdzcil,
      Cpar,   Vp,    Cfs,   Cdz,   Cdzil,
      Cfsil,  Cfs,   Vp,    Cfs,   Cfsil,
      Cdzil,  Cdz,   Cfs,   Vp,    Cpar,
      Cdzcil, Cdzil, Cfsil, Cpar,  Vp
    ),
    5, 5,
    byrow = T
  )
  
  # phenotypic means
  me <- rep(0, 5)
  
  # true parameter
  true <- data.frame(
    est = c(SA1, SC1, ST1, SD1, SI1, SEW1, SEB1, SE1),
    component = c("A", "C", "T", "D", "N",  "EW", "EB", "E")
  ) # note the difference due to the increase in phenotype variance given q as a consequence of AM
  
  # iterate simulation and estimate parameters
  est_ctd  <- c()
  sel_vars <- c("stw1", "tw1", "s1", "tw2", "stw2")
  sel_vars_tw <- c("tw1", "tw2")
  
  ## exact data simulation ####
  zyg <- c(rep(1, Nmz), rep(2, Ndz))
  dat <- matrix(0, N, Nmem1)
  dat[, 1] <- zyg
  dat[zyg == 1, 2:Nmem1] <- mvrnorm(Nmz, mu = me, Sigma = Smz / cor_Nmz, emp = exact) # note here is sample parameters
  dat[zyg == 2, 2:Nmem1] <- mvrnorm(Ndz, mu = me, Sigma = Sdz / cor_Ndz, emp = exact) # note here is sample parameters
  colnames(dat) <- c("zyg", sel_vars)
  dat <- as.data.frame(dat)
  
  # data files
  mz <- dat[dat$zyg == 1, c("zyg", sel_vars_tw)][, -1]
  dz <- dat[dat$zyg == 2, c("zyg", sel_vars_tw)][, -1]
  mzfull <- dat[dat$zyg == 1, c("zyg", sel_vars)][, -1]
  dzfull <- dat[dat$zyg == 2, c("zyg", sel_vars)][, -1]
  
  # extract correlations
  rmz <- cor(mz)["tw1", "tw2"]
  rdz <- cor(dz)["tw1", "tw2"]
  rsib <- cor(mzfull)["tw1", "s1"]
  rmzil <- cor(mzfull)["tw1", "stw2"]
  rdzil <- cor(dzfull)["tw1", "stw2"]
  rsibil <- cor(mzfull)["s1", "stw2"]
  rmzcil <- cor(mzfull)["stw1", "stw2"]
  rdzcil <- cor(dzfull)["stw1", "stw2"]
  ram <- cor(mzfull)["tw1", "stw1"]
  
  return(
    if(raw == T){
      
      list(
        dat   = dat
      )
    } else {
      
      list(value = list(Vp = Vp, Vs = Vs, rmz = rmz, rsib = rsib, rmzil = rmzil, rdzil = rdzil, rsibil = rsibil, rmzcil = rmzcil, rdzcil = rdzcil, ram = ram))
    }    
  )
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #    
}

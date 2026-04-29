# ADE or ACE function
# Written by G.B. (Giaco)
# This script is based on H.H.M. CTD script https://hermine-maes.squarespace.com/

require(OpenMx)
acde.fun <- function(
    
    mzData=NULL,
    dzData=NULL,
    selVars = c("T1_p","T2_p"),
    name = "ADCE",
    ...,
    
    # Phenotype mean
    M =  1,
    
    # If any mean effect is expected change label accordingly / default is one mean
    # Labels for mean MZ 
    mean.MZ = c("mean", "mean"),
    # Labels for mean DZ
    mean.DZ = c("mean", "mean"),
    
    # Phenotype variances
    # by default fit an AE model
    ap.free = T, # Freely estimate a? 
    cp.free = F, # Freely estimate c? 
    dp.free = F, # Freely estimate d? 
    
    ap.val  = .30, # Starting value for a
    cp.val  = .00, # Starting value for c
    dp.val  = .00, # Starting value for d
    ep.val  = .70  # Starting value for e
    
){
  # ----------------------------------------------------------------------------------------------------------------------
  # Create Algebra for expected Mean Matrices
  meanGMZ <- mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=M, labels=c("meanMZ1", "meanMZ2"), name="meanGMZ" )
  meanGDZ <- mxMatrix( type="Full", nrow=1, ncol=2, free=TRUE, values=M, labels=c("meanDZ1", "meanDZ2"), name="meanGDZ" )
  
  # Create Matrices for Variance Components
  covA <- mxMatrix( type="Symm", nrow=1, ncol=1, free=ap.free, values=ap.val, label="VA11", name="VA" )
  covC <- mxMatrix( type="Symm", nrow=1, ncol=1, free=cp.free, values=cp.val, label="VC11", name="VC" )
  covD <- mxMatrix( type="Symm", nrow=1, ncol=1, free=dp.free, values=dp.val, label="VD11", name="VD" )
  covE <- mxMatrix( type="Symm", nrow=1, ncol=1, free=TRUE,    values=ep.val, label="VE11", name="VE" )
  
  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covP  <- mxAlgebra( expression = VA + VD + VC + VE, name="V" )
  covMZ <- mxAlgebra( expression = VA + VD + VC, name="cMZ" )
  covDZ <- mxAlgebra( expression = 0.5%x%VA + 0.25%x%VD + VC, name="cDZ" )
  expCovMZ <- mxAlgebra( expression = rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
  expCovDZ <- mxAlgebra( expression = rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
  
  # Create Data Objects for Multiple Groups
  dataMZ <- mxData( observed=mzData, type="raw" )
  dataDZ <- mxData( observed=dzData, type="raw" )
  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal( covariance="expCovMZ", means="meanGMZ", dimnames=selVars )
  expDZ <- mxExpectationNormal( covariance="expCovDZ", means="meanGDZ", dimnames=selVars )
  funML <- mxFitFunctionML()
  # Create Model Objects for Multiple Groups
  pars <- list( covA, covD, covC, covE, covP )
  modelMZ <- mxModel( pars, meanGMZ, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
  modelDZ <- mxModel( pars, meanGDZ, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
  multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
  # Create Algebra for Unstandardized and Standardized Variance Components
  rowUS <- c('US')
  colUS <- rep(c('varP','varP_A','varP_D','varP_C','varP_E'),each=1)
  estUS <- mxAlgebra( expression=cbind(V,VA/V,VD/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
  # Create Confidence Interval Objects
  ciADE <- mxCI( "US[1,1:5]" )
  # Build Model with Confidence Intervals
  modelACDE <- mxModel( "oneADEvc", pars, modelMZ, modelDZ, multi, estUS, ciADE )
  
}

#
# Code supporting
#	Metacommunity dynamics and the spurious detection of species associations in co-occurrence analyses
#	Vincent Calcagno, Nik J Cunniffe and Frederic M Hamelin
#
# The code will simulate patch occupancy matrices according to a generalised metacommunity model
#	in which species do not interact, all patches are identical and in which homogeneous dispersal is assumed,
#  	but in which patch disturbance occurs at a (user-defined) rate depending on the patch age.
#
# The code also does various calculations - as outlined in the paper - allowing predictions to be made of which
#	pairs of species will spuriously be assumed to be associated when using common tests of species assocation
# 	as well as wrapping up the test using the permutation algorithms Sim2 and Sim9 introduced by Gotelli (and co-authors)
#
# The associated R Markdown file shows how the code can be used in context, providing a repository generating all figures in the paper
#

#
# For calculating the Ulrich and Gotelli metrics
#
library(EcoSimR)

#
# For accurate numerical integration
#
library(cubature)

######################
# Extenal functions  #
######################

#
# Return an empty list of model information
#
# Parameters largestMaxAge and intTol control how integration is done internally
#
createModel <- function(largestMaxAge = 1e4,   # (Large) numerical value that is used if patch max age is not set and so if ever need to integrate to "infinity"
                        intTol        = 1e-5,  # Relative tolerance for numerical integration
                        maxIts        = 250,   # Maximum number of iterations in the procedure to calculate model equilibria
                        sqChange      = 1e-25) # Termination condition the procedure to calculate model equilibria
{
  newModel <- list(
                    internalConstants = list(integrationMaxAge = largestMaxAge, integrationTol = intTol, eqmMaxIterations = maxIts, eqmSqChange = sqChange),
                    modelFlags        = list(setParameters = FALSE, setPatchDeath = FALSE, doneCalculations = FALSE),
                    modelParameters   = list(numSpecies = NA, vM = NA, vC = NA, vE = NA),
                    modelPatchDeath   = list(deathRateFunc = NA, integratedDeathRateFunc = NA, maxAge = NA, type = NA, r = NA, b = NA, vPts = NA, vVals = NA, nPts = NA,scaleFact = NA),
                    modelCalculations = list(vA = NA, vR = NA, vP = NA, mCovariance = NA, vVariance = NA, averageFastness = NA)
                  )
  return(newModel)
}

#
# Sets model parameters
#
# Parameters n,M,C and E set number of species and rates in the metapopulation model (apart from death rate)
#
initParameters <- function(thisModel, 	  # Model structure, as previously returned by createModel()
                               n,         # Number of species
                               M,         # Vector of migration rates
                               C,         # Vector of colonisation rates
                               E)         # Vector of extinction rates
{
  # setting or altering parameters will invalidate any calculations
  thisModel$modelFlags$doneCalculations <- FALSE
  thisModel$modelParameters$numSpecies <- n
  thisModel$modelParameters$vM <- M
  thisModel$modelParameters$vC <- C
  thisModel$modelParameters$vE <- E
  thisModel$modelFlags$setParameters <- TRUE
  return(thisModel)
}

#
# Sets up the general model for death of patches
#
initPatchDeath <- function(thisModel,       	 # Model structure, as previously returned by createModel()
                               deathRateFunc = internalDefaultDeathRateFunc,
                               integratedDeathRateFunc = internalDefaultIntegratedDeathRateFunc, # pass NULL for numerical integration of deathRateFunc
                               maxAge = Inf,     # Maximum ages of patches
                               type = 1,         # (If using defaultDeathRateFunc, which type to use)
                               r = 1,            # (If using defaultDeathRateFunc, which value of r to use)
                               b = NA,           # (If using defaultDeathRateFunc, which value of b to use)
                               numPoints = 1000, # (If using numerical integration of deathRateFunc) minimum number of cached points
                               expThresh = 0.01) # (If using numerical integration of deathRateFunc) maximum by which exp(-deathRateFunc) can change between two successive samples
{
  # setting or altering death rates will invalidate any calculations
  thisModel$modelFlags$doneCalculations <- FALSE
  thisModel$modelPatchDeath$maxAge <- maxAge
  thisModel$modelPatchDeath$deathRateFunc <- deathRateFunc
  if(!is.null(integratedDeathRateFunc))
  {
    thisModel$modelPatchDeath$integratedDeathRateFunc <- integratedDeathRateFunc
  }else{
    thisModel <- internalSetupNumericalIntegratedDeathRate(thisModel,numPoints,expThresh)
    thisModel$modelPatchDeath$integratedDeathRateFunc <- internalNumericallyIntegratedDeathRateFunc
  }
  thisModel$modelPatchDeath$type <- type
  thisModel$modelPatchDeath$r <- r
  thisModel$modelPatchDeath$b <- b
  f <- function(x,thisModel) exp(-thisModel$modelPatchDeath$integratedDeathRateFunc(x,thisModel))
  intMaxAge <- min(thisModel$internalConstants$integrationMaxAge,thisModel$modelPatchDeath$maxAge)
  v <- hcubature(f,
                 lower=0,
                 upper=intMaxAge,
                 tol=thisModel$internalConstants$integrationTol,
                 thisModel=thisModel)
  thisModel$modelPatchDeath$scaleFact <- v$integral
  thisModel$modelFlags$setPatchDeath <- TRUE
  return(thisModel)
}

doCalculations <- function(thisModel)
{
  if(thisModel$modelFlags$setParameters & thisModel$modelFlags$setPatchDeath)
  {
    thisModel$modelCalculations$vP <- internalFindEqm(thisModel)

    thisModel$modelCalculations$vR <- thisModel$modelParameters$vM + thisModel$modelParameters$vC * thisModel$modelCalculations$vP + thisModel$modelParameters$vE
    thisModel$modelCalculations$vA <- (thisModel$modelParameters$vM + thisModel$modelParameters$vC * thisModel$modelCalculations$vP) / thisModel$modelCalculations$vR

    thisModel$modelCalculations$mCovariance <- internalFindCovarianceMatrix(thisModel)
    thisModel$modelCalculations$vVariance <- diag(thisModel$modelCalculations$mCovariance)
    thisModel$modelCalculations$averageFastness <- internalAverageFastness(thisModel)

    thisModel$modelFlags$doneCalculations <- TRUE
  }
  return(thisModel)
}

#
# Take a sample from the community matrix according to a given model
#
sampleCommunityMatrix <- function(thisModel,nPatches)
{
  retVal <- list(cMatrix = NA, patchAges = NA)
  if(thisModel$modelFlags$doneCalculations)
  {
    # sample patch ages
    vAge <- internalSamplePatchAge(nPatches,thisModel)
    retVal$patchAges <- vAge

    # all patches start off empty of all species
    retVal$cMatrix <- matrix(0,nPatches,thisModel$modelParameters$numSpecies)
    for(i in 1:nPatches)
    {
      for(j in 1:thisModel$modelParameters$numSpecies)
      {
        thisP <- thisModel$modelCalculations$vA[j] * (1 - exp(-thisModel$modelCalculations$vR[j]*vAge[i]))
        # if appropriate, fill in species j in patch i
        if(runif(1)<thisP)
        {
          retVal$cMatrix[i,j]<-1
        }
      }
    }
  }
  return(retVal)
}

#
# Initialise the pairwiseC information by doing the calculation on the data
#
initPairwise <- function(commMatrix)
{
  nSpp <- dim(commMatrix)[2]
  newPairwise <- list(
    data           = list(commMatrix = commMatrix, nSpp = nSpp),
    speciesPairs   = list(nPairs = NA, allSppOne = NA, allSppTwo = NA),
    dataPairwiseC  = list(C = NA),
    simPairwiseC   = list(pairwiseScores = NA),
    output         = list(pHigh = NA, pLow = NA, stEffSize = NA, pHighAdj = NA, pLowAdj = NA)
  )
  newPairwise$speciesPairs$nPairs    <- nSpp*(nSpp-1)/2
  newPairwise$speciesPairs$allSppOne <- rep(NA,newPairwise$speciesPairs$nPairs)
  newPairwise$speciesPairs$allSppTwo <- rep(NA,newPairwise$speciesPairs$nPairs)
  thisIDX <- 1
  for(spOne in 1:nSpp)
  {
    spTwo <- spOne + 1
    while(spTwo <= nSpp)
    {
      newPairwise$speciesPairs$allSppOne[thisIDX] <- spOne
      newPairwise$speciesPairs$allSppTwo[thisIDX] <- spTwo
      thisIDX <- thisIDX + 1
      spTwo <- spTwo + 1
    }
  }
  newPairwise$dataPairwiseC$C <- internalPairwiseC(newPairwise$data$commMatrix,newPairwise$speciesPairs$allSppOne,newPairwise$speciesPairs$allSppTwo)
  return(newPairwise)
}

#
# Finalise the pairwiseC information by doing the simulations and the calculations on the results
#
simPairwise <- function(pairwiseInfo,simReps,burnInReps,doSim9)
{
  #
  # do the reshuffling
  #
  speciesData <- t(pairwiseInfo$data$commMatrix)

  metricF <- get("c_score")
  Obs <- metricF(speciesData)
  msim <- speciesData[rowSums(speciesData) > 0, ]
  burn.in.metric <- vector(mode = "numeric", length = burnInReps)
  simulated.metric <- vector(mode = "numeric", length = simReps)
  for (i in 1:burnInReps) {
    if(doSim9){
      msim <- sim9_single(msim)
    }else{
      msim <- sim2(msim)
    }
    burn.in.metric[i] <- metricF(msim)
  }
  pairwiseInfo$simPairwiseC$pairwiseScores <- matrix(NA,nrow=pairwiseInfo$speciesPairs$nPairs,ncol=simReps)
  for (i in 1:simReps) {
    if(doSim9){
      msim <- sim9_single(msim)
    }else{
      msim <- sim2(msim)
    }
    simulated.metric[i] <- metricF(msim)
    v <- internalPairwiseC(t(msim),
                           pairwiseInfo$speciesPairs$allSppOne,
                           pairwiseInfo$speciesPairs$allSppTwo)
    pairwiseInfo$simPairwiseC$pairwiseScores[,i] <- v
  }
  #
  # Do pairwise comparisons of simulation results against the data
  #
  pairwiseInfo$output$pHigh <- numeric(pairwiseInfo$speciesPairs$nPairs)
  pairwiseInfo$output$pLow  <- numeric(pairwiseInfo$speciesPairs$nPairs)
  pairwiseInfo$output$stEffSize  <- numeric(pairwiseInfo$speciesPairs$nPairs)
  for(i in 1:pairwiseInfo$speciesPairs$nPairs)
  {
    # to do statistical testing, there are two p-values we might care about

    # pHigh is p-value relevant to testing whether the value of c_ij is
    # in the upper tail and so is higher than expected (i.e. segregation)

    # pLow is p-value relevant to testing whether the value of c_ij is
    # in the lower tail and so is lower than expected (i.e. aggregation)
    pairwiseInfo$output$pHigh[i] <- sum(pairwiseInfo$simPairwiseC$pairwiseScores[i,]>pairwiseInfo$dataPairwiseC$C[i])/length(pairwiseInfo$simPairwiseC$pairwiseScores[i,])
    pairwiseInfo$output$pLow[i] <- sum(pairwiseInfo$simPairwiseC$pairwiseScores[i,]<pairwiseInfo$dataPairwiseC$C[i])/length(pairwiseInfo$simPairwiseC$pairwiseScores[i,])

    # another quantity that might be of interest is the standardised effect size
    # https://esajournals.onlinelibrary.wiley.com/doi/epdf/10.1002/ecs2.2797
    pairwiseInfo$output$stEffSize[i] <- (pairwiseInfo$dataPairwiseC$C[i] - mean(pairwiseInfo$simPairwiseC$pairwiseScores[i,]))/sd(pairwiseInfo$simPairwiseC$pairwiseScores[i,])
  }

  #
  # adjust p-values to account for doing so many comparisons
  #
  pairwiseInfo$output$pHighAdj <- p.adjust(pairwiseInfo$output$pHigh, method = "fdr")
  pairwiseInfo$output$pLowAdj <- p.adjust(pairwiseInfo$output$pLow, method = "fdr")

  return(pairwiseInfo)
}

######################
# Internal functions #
######################

internalPairwiseC <- function(commMatrix,allSppOne,allSppTwo)
{
  n <- length(allSppOne)
  vC <- numeric(n)
  for(i in 1:n)
  {
    rOne <- sum(commMatrix[,allSppOne[i]])
    rTwo <- sum(commMatrix[,allSppTwo[i]])
    S <- sum((commMatrix[,allSppOne[i]]+commMatrix[,allSppTwo[i]])==2)
    thisC <- (rOne - S)*(rTwo - S)
    vC[i] <- thisC
  }
  return(vC)
}

internalSetupNumericalIntegratedDeathRate <- function(thisModel,numPoints,eThresh)
{
  maxAge <- min(thisModel$internalConstants$integrationMaxAge,thisModel$modelPatchDeath$maxAge)
  # first do a random sample of numPoints in between 0 and maxAge
  thisModel$modelPatchDeath$vPts <- sort(c(0,runif(numPoints,0,maxAge),maxAge,Inf))  # artificially inflate number of points
                                                                                     # to ensure later look up works smoothly
  thisModel$modelPatchDeath$nPts <- numPoints+3
  thisModel$modelPatchDeath$vVals <- numeric(thisModel$modelPatchDeath$nPts)
  thisModel$modelPatchDeath$vVals[1] <- 0
  for(i in 2:(thisModel$modelPatchDeath$nPts-1))
  {
    xStart <- thisModel$modelPatchDeath$vPts[i-1]
    xEnd <- thisModel$modelPatchDeath$vPts[i]
    v <- hcubature(thisModel$modelPatchDeath$deathRateFunc,
                   lower=xStart,
                   upper=xEnd,
                   tol=thisModel$internalConstants$integrationTol,
                   thisModel=thisModel)
    thisModel$modelPatchDeath$vVals[i] <- thisModel$modelPatchDeath$vVals[i-1] + v$integral
  }
  thisModel$modelPatchDeath$vVals[thisModel$modelPatchDeath$nPts] <- thisModel$modelPatchDeath$vVals[thisModel$modelPatchDeath$nPts-1]
  #
  # add additional sampling points to ensure exp(-deathRateFunc) never changes too quickly between points
  #
  vNewPts <- numeric(thisModel$modelPatchDeath$nPts)
  vNewVals <- numeric(thisModel$modelPatchDeath$nPts)
  vNewPts[1] <- thisModel$modelPatchDeath$vPts[1]
  vNewVals[1] <- thisModel$modelPatchDeath$vVals[1]
  lastX <- vNewPts[1]
  lastY <- vNewVals[1]
  lastE <- exp(-lastY)
  fromIDX <- 2
  nextIDX <- 1
  nPts <- thisModel$modelPatchDeath$nPts
  while(fromIDX <= thisModel$modelPatchDeath$nPts)
  {
    trialX <- thisModel$modelPatchDeath$vPts[fromIDX]
    trialY <- thisModel$modelPatchDeath$vVals[fromIDX]
    trialE <- exp(-trialY)
    tookFrom <- TRUE
    while((lastE-trialE)>=eThresh)
    {
      trialX <- lastX + (trialX - lastX)/2.0
      v <- hcubature(thisModel$modelPatchDeath$deathRateFunc,
                     lower=lastX,
                     upper=trialX,
                     tol=thisModel$internalConstants$integrationTol,
                     thisModel=thisModel)
      trialY <- lastY + v$integral
      trialE <- exp(-trialY)
      tookFrom <- FALSE
    }
    nextIDX <- nextIDX + 1
    vNewPts[nextIDX] <- trialX
    vNewVals[nextIDX] <- trialY
    lastX <- trialX
    lastY <- trialY
    lastE <- trialE
    if(tookFrom)
    {
      fromIDX <- fromIDX + 1
    }
  }
  thisModel$modelPatchDeath$vPts <- vNewPts
  thisModel$modelPatchDeath$vVals <- vNewVals
  thisModel$modelPatchDeath$nPts <- nextIDX
  return(thisModel)
}

internalNumericallyIntegratedDeathRateFunc <- function(age,thisModel)
{
  # make sure the function can be called in a vectorised context
  n <- length(age)
  # default to return 0 unless age is valid
  retVal <- rep(0,n)
  lastAge <- NA
  for(i in 1:n)
  {
    if(age[i] < 0)
    {
      retVal[i] <- 0
    }else
    {
      if(age[i] <= thisModel$modelPatchDeath$maxAge)
      {
        thisAge <- age[i]
      }else{
        thisAge <- thisModel$modelPatchDeath$maxAge
      }
      belowIDX <- max(which(thisModel$modelPatchDeath$vPts<=thisAge))
      if(length(which(thisModel$modelPatchDeath$vPts>thisAge))>0)
      {
        aboveIDX <- min(which(thisModel$modelPatchDeath$vPts>thisAge))
        belowX <- thisModel$modelPatchDeath$vPts[belowIDX]
        aboveX <- thisModel$modelPatchDeath$vPts[aboveIDX]
        belowY <- thisModel$modelPatchDeath$vVals[belowIDX]
        aboveY <- thisModel$modelPatchDeath$vVals[aboveIDX]
        if(aboveIDX < thisModel$modelPatchDeath$nPts)
        {
          thisVal <- belowY + ((thisAge - belowX)/(aboveX - belowX))*(aboveY - belowY)
        }else{
          # simple interpolation not work correctly if upper limit is Inf, so handle separately
          thisVal <- aboveY
        }
      }else{
        thisVal <- thisModel$modelPatchDeath$vVals[thisModel$modelPatchDeath$nPts]
      }
      retVal[i] <- thisVal
    }
  }
  return(retVal)
}

#
# Prepackaged death rates
#
internalDefaultDeathRateFunc <- function(age,thisModel)
{
  # make sure the function can be called in a vectorised context
  n <- length(age)
  # default to return 0 unless age is valid
  retVal <- rep(0,n)
  for(i in 1:n)
  {
    if(age[i] >= 0 & age[i] <= thisModel$modelPatchDeath$maxAge)
    {
      if(thisModel$modelPatchDeath$type==1)
      {
        # constant
        retVal[i] <- thisModel$modelPatchDeath$r
      }
      if(thisModel$modelPatchDeath$type==2)
      {
        # exponentially decreasing
        retVal[i] <- thisModel$modelPatchDeath$r*exp(-age[i]/thisModel$modelPatchDeath$b)
      }
    }
  }
  return(retVal)
}

#
# Prepackaged integrated death rates
#
internalDefaultIntegratedDeathRateFunc <- function(age,thisModel)
{
  # make sure the function can be called in a vectorised context
  n <- length(age)
  retVal <- rep(NA,n)
  for(i in 1:n)
  {
    if(age[i] < 0)
    {
      retVal[i] <- 0
    }else
    {
      if(age[i] <= thisModel$modelPatchDeath$maxAge)
      {
        thisAge <- age[i]
      }else{
        thisAge <- thisModel$modelPatchDeath$maxAge
      }
      if(thisModel$modelPatchDeath$type==1)
      {
        # constant
        retVal[i] <- thisModel$modelPatchDeath$r * thisAge
      }
      if(thisModel$modelPatchDeath$type==2)
      {
        # exponentially decreasing
        retVal[i] <- thisModel$modelPatchDeath$r*thisModel$modelPatchDeath$b*(1 - exp(-thisAge/thisModel$modelPatchDeath$b))
      }
    }
  }
  return(retVal)
}

#
# Probability density function for patch ages
#   (it is called $p_{(\cdot,x}$ in the paper)
#
# Note it is scaled to be a probability density function
#
internalPdfPatchAge <- function(age,thisModel)
{
  retVal <- exp(-thisModel$modelPatchDeath$integratedDeathRateFunc(age,thisModel))
  retVal <- retVal / thisModel$modelPatchDeath$scaleFact
  return(retVal)
}

internalPdfPatchAgeIntegrated <- function(age,thisModel)
{
  n <- length(age)
  retVal <- numeric(n)
  for(i in 1:n)
  {
    v <- hcubature(internalPdfPatchAge,
                 lower=0,
                 upper=age[i],
                 tol=thisModel$internalConstants$integrationTol,
                 thisModel=thisModel)
    retVal[i] <- v$integral
  }
  return(retVal)
}

#
# Probability density function for patch being of a certain age and containing species controlled by trialA and trialR
#
internalF <- function(age,trialA,trialR,thisModel)
{
  p_ix <- trialA * (1 - exp(-trialR*age))
  p_dotx <- internalPdfPatchAge(age,thisModel)
  return(p_ix*p_dotx)
}

#
# Iterative routine to find the equilibrium patch densities
#
internalFindEqm <- function(thisModel)
{
  maxAge <- min(thisModel$internalConstants$integrationMaxAge,thisModel$modelPatchDeath$maxAge)
  vTrialP <- rep(1/2,thisModel$modelParameters$numSpecies)
  k <- 1
  sqDiff <- 1
  while(k <= thisModel$internalConstants$eqmMaxIterations & sqDiff >= thisModel$internalConstants$eqmSqChange)
  {
    vTrialR <- thisModel$modelParameters$vM + thisModel$modelParameters$vC * vTrialP + thisModel$modelParameters$vE
    vTrialA <- (thisModel$modelParameters$vM + thisModel$modelParameters$vC * vTrialP) / vTrialR
    for(i in 1:thisModel$modelParameters$numSpecies)
    {
      v <- hcubature(internalF,
                     lower=0,
                     upper=maxAge,
                     trialA=vTrialA[i],
                     trialR=vTrialR[i],
                     thisModel=thisModel,
                     tol=thisModel$internalConstants$integrationTol)
      vTrialP[i] <- v$integral
    }
    if(k > 1)
    {
      sqDiff <- sum((vPLast-vTrialP)^2)
    }
    vPLast <- vTrialP
    k <- k + 1
  }
  return(vTrialP)
}

#
# Probability density function for patch of age x containing sppOne and sppTwo simultaneously
#
internalH <- function(age,sppOne,sppTwo,thisModel)
{
  pi_ix <- thisModel$modelCalculations$vA[sppOne]*(1 - exp(-thisModel$modelCalculations$vR[sppOne]*age))/thisModel$modelCalculations$vP[sppOne]
  pi_jx <- thisModel$modelCalculations$vA[sppTwo]*(1 - exp(-thisModel$modelCalculations$vR[sppTwo]*age))/thisModel$modelCalculations$vP[sppTwo]
  p_dotx <- internalPdfPatchAge(age,thisModel)
  return(pi_ix*pi_jx*p_dotx)
}

#
# Probability density function for patch of age x containing thisSpp
#
internalG <- function(age,thisSpp,thisModel)
{
  pi_ix <- thisModel$modelCalculations$vA[thisSpp]*(1 - exp(-thisModel$modelCalculations$vR[thisSpp]*age))/thisModel$modelCalculations$vP[thisSpp]
  p_dotx <- internalPdfPatchAge(age,thisModel)
  return(pi_ix*p_dotx)
}

#
# Calculates covariance between species (integrating over all patch ages)
#
internalFindCovarianceMatrix <- function(thisModel)
{
  cMatrix <- matrix(0,nrow=thisModel$modelParameters$numSpecies,ncol=thisModel$modelParameters$numSpecies)
  maxAge <- min(thisModel$internalConstants$integrationMaxAge,thisModel$modelPatchDeath$maxAge)
  for(i in 1:thisModel$modelParameters$numSpecies)
  {
    for(j in i:thisModel$modelParameters$numSpecies)
    {
      v1 <- hcubature(internalH,
                      lower=0,
                      upper=maxAge,
                      sppOne=i,
                      sppTwo=j,
                      thisModel=thisModel,
                      tol=thisModel$internalConstants$integrationTol)
      v2 <- hcubature(internalG,
                      lower=0,
                      upper=maxAge,
                      thisSpp=i,
                      thisModel=thisModel,
                      tol=thisModel$internalConstants$integrationTol)
      cMatrix[i,j] <- v1$integral - v2$integral*v2$integral # note should always be 1
      cMatrix[j,i] <- cMatrix[i,j]
    }
  }
  return(cMatrix)
}

#
# Average profile of species, weighted by equilibrium abundance
#
internalW <- function(age,thisModel)
{
  vW <- thisModel$modelCalculations$vP/sum(thisModel$modelCalculations$vP)
  n <- length(age)
  retVal <- rep(NA,n)
  for(k in 1:n)
  {
    retVal[k] <- 0
    for(i in 1:thisModel$modelParameters$numSpecies)
    {
      thisW <- vW[i]*thisModel$modelCalculations$vA[i]*(1 - exp(-thisModel$modelCalculations$vR[i]*age[k]))/thisModel$modelCalculations$vP[i]
      retVal[k] <- retVal[k] + thisW
    }
  }
  return(retVal)
}

#
# Probability density function for patch of age x containing "average" species
#
internalJ <- function(age,thisModel)#x,nSpp,vA,vR,vP,myParams,deathRateFunc,integratedDeathRateFunc)
{
  pi_x <- internalW(age,thisModel)
  p_dotx <- internalPdfPatchAge(age,thisModel)
  return(pi_x*p_dotx)
}

#
#  Probability density function for squared density in patch of age x containing "average" species
#
internalK <- function(age,thisModel)
{
  pi_x <- internalW(age,thisModel)
  p_dotx <- internalPdfPatchAge(age,thisModel)
  return(pi_x*pi_x*p_dotx)
}

#
# Calculate the average "fastness" (i.e. the variance of the average species)
#   (this is $\tilde \pi_x$ in the notes)
#
internalAverageFastness <- function(thisModel)
{
  maxAge <- min(thisModel$internalConstants$integrationMaxAge,thisModel$modelPatchDeath$maxAge)

  v1 <- hcubature(internalK,
                  lower=0,
                  upper=maxAge,
                  thisModel=thisModel,
                  tol=thisModel$internalConstants$integrationTol)

  v2 <- hcubature(internalJ,
                  lower=0,
                  upper=maxAge,
                  thisModel=thisModel,
                  tol=thisModel$internalConstants$integrationTol)

  c <- v1$integral - v2$integral*v2$integral # note  should always be 1
  return(c)
}

#
# Sample patch ages according to the death rate
#
internalSamplePatchAge <- function(nPatches,thisModel)
{
  # do patches have maximum age?
  maxAge <- min(thisModel$internalConstants$integrationMaxAge,thisModel$modelPatchDeath$maxAge)
  allSamples <- numeric(nPatches)
  #
  # sample out of patch age density function
  #
  for(i in 1:nPatches)
  {
    d <- runif(1)
    #
    # Now want to solve integrateDeathRate(x) = -log(1-d) for x
    #
    f <- function(x) internalPdfPatchAgeIntegrated(x,thisModel) - d
    allSamples[i] <- uniroot(f,c(0,maxAge))$root
  }
  return(allSamples)
}


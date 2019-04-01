############################################################################
## Ordination code for final paper in BIOL 419: Modern Statistical Models ##
## Author: Vaughn Shirey, April 1st, 2019                                 ##
############################################################################

## Import libraries ##
library(ggfortify)

## Set seed ##
set.seed(04012019)

## Establish constants ##
numSpecies <- 3
numSites <- 3

## Establish environmental variable distributions for sampling ##


## Establish simulation function ##
sim <- function(numInd = 20, numTrt = 5, trtVar = "normal", missDat = FALSE){
  #################################################################################################################
  ## Simulates morphometric trait data from a population of three (3) butterflies along an elevational gradient. ##
  ##                                                                                                             ##
  ## numInd = number of individuals to sample per species.                                                       ##
  ## numTrt = number of traits to sample per individual.                                                         ##
  ## trtVar = variance within a particular trait ("low", "normal", "high", "mixed").                             ##
  ## missDat = should missing data be injected into the final dataset.                                           ##
  #################################################################################################################
  
  ## Generate environmental variables per elevation ##
  loElevTemp <- rep(rnorm(1, mean=20, sd=1), numInd)
  midElevTemp <- rep(rnorm(1, mean=15, sd=1), numInd)
  hiElevTemp <- rep(rnorm(1, mean=5, sd=1), numInd)
    
  loElevWind <- rep(abs(rnorm(1, mean=5, sd=2)), numInd)
  midElevWind <- rep(abs(rnorm(1, mean=10, sd=2)), numInd)
  hiElevWind <- rep(abs(rnorm(1, mean=10, sd=2)), numInd)
  
  siteCode1 <- rep(1, numInd)
  siteCode2 <- rep(2, numInd)
  siteCode3 <- rep(3, numInd)
  
  site <- append(append(siteCode1, siteCode2), siteCode3)
  
  temp <- append(append(loElevTemp, midElevTemp), hiElevTemp)
  wind <- append(append(loElevWind, midElevWind), hiElevWind)
  
  ## Establish environmental matrix ##
  envMatrix <- matrix(nrow = numInd*numSites, ncol = 3)
  envMatrix[,1] <- site
  envMatrix[,2] <- temp
  envMatrix[,3] <- wind
  
  ## Generate trait variables per individual per site ##
  tempMatrix <- matrix(nrow = numInd*numSites, ncol=numTrt)
  
  if(trtVar == "normal"){
    for(i in 1:numTrt){
      loTrt <- ((loElevTemp[1]+loElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=1))
      midTrt <- ((midElevTemp[1]+midElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=1))
      hiTrt <- ((hiElevTemp[1]+hiElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=1))
      
      tempMatrix[,i] <- append(append(loTrt, midTrt), hiTrt)
    }
  }
  
  if(trtVar == "low"){
    for(i in 1:numTrt){
      loTrt <- ((loElevTemp[1]+loElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=0.5))
      midTrt <- ((midElevTemp[1]+midElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=0.5))
      hiTrt <- ((hiElevTemp[1]+hiElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=0.5))
      
      tempMatrix[,i] <- append(append(loTrt, midTrt), hiTrt)
    }
  }
  
  if(trtVar == "high"){
    for(i in 1:numTrt){
      loTrt <- ((loElevTemp[1]+loElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=1.5))
      midTrt <- ((midElevTemp[1]+midElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=1.5))
      hiTrt <- ((hiElevTemp[1]+hiElevWind[1])*rnorm(1, mean=1, sd=0.25))*abs(rnorm(numInd, mean=0, sd=1.5))
      
      tempMatrix[,i] <- append(append(loTrt, midTrt), hiTrt)
    }
  }
  
  ## Establish trait matrix ##
  traitMatrix <- matrix(nrow = numInd*numSites, ncol=numTrt+1)
  traitMatrix[,1] <- site
  traitMatrix[,1:numTrt+1] <- tempMatrix
    
  return(list(envMatrix, traitMatrix))
}

data <- sim()
df.traits <- as.data.frame(data[[2]])
df.traits$V1 <- as.factor(df.traits$V1)
df.env <- as.data.frame(data[[1]])

## Principal Components Analysis ##
time.init <- Sys.time()
autoplot(prcomp(df.traits[,2:ncol(df.traits)]), data=df.traits, colour='V1', loadings = TRUE, loadings.label = TRUE)
pca.Time <- Sys.time() - time.init

## Detrended Correspondance Analysis ##

## Canonical Correspondance Analysis ##

## Non-metric Multidimensional Scaling ##


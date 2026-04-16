library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <- 15 #Number of cores
set.seed(150)
cluster = makeCluster(iH)
vSeeds = sample(1:1e5, iH, replace = FALSE)
clusterExport(cluster, c("vSeeds"))
clusterEvalQ(cluster, {
  library(RandAlphaTest)
})
for (h in 1:iH) {
  clusterSetRNGStream(cl = cluster, vSeeds[h])
}
iM = 1000 # Number of MC samples


vN = c(100, 200, 500) # Cross-sectional sizes
vT = c(100, 200, 300, 500) # Temporal sizes

# Load simulated data. Change according to your working directory
sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Power/"
sType = "power"
load(paste(sPathData, sType, "LatentFactor_PhiNu040.Rdata", sep = ""))

dAlpha = 0.05 # Level of the one-shot tests
dNu <- 5 # Power for psi
dQ = 1/4 # Power for f(B) in the modified threshold

# Pre-allocation matrices
mRes1 = matrix(0.0, ncol = length(vT), nrow = length(vN)) # For LIL-based threshold 
mRes2 = matrix(0.0, ncol = length(vT), nrow = length(vN)) # For f(B) based threshold

for(n  in 1:length(vN)){
  iN = vN[n]
  iB = (floor(log(iN)^2)) # Number of trials in the de-randomized procedure

  mZ = matrix(0.0, nrow = iM, ncol =iB)

  lAB = get_AB(iN) # normalizing sequences for EVT asymptotics
  dC = get_criticalValues_ultimate(lAB$A, lAB$B, dAlpha) # Gumbel-based critical value
  lFoo = lData[[n]]

  dCritic1 = (1.0-dAlpha) - sqrt(dAlpha*(1.0-dAlpha))*sqrt(2.0*log(log(iB))/iB) # LIL-based threshold
  dCritic2 = (1.0-dAlpha) - iB^(-dQ)  # f(B) based threshold

  for(j in 1:length(vT)){
    iT = vT[j]
    dTeff <- iT^(1/dNu)

    for(m in 1:iM){
      set.seed(m)
      
      lFoo1 = list(Y = lFoo[[m]]$Y[, 1:iT], X = lFoo[[m]]$X[, 1:iT])
      clusterExport(cl=cluster, c("dNu", "dTeff"))

      # Test statistics for the B trials over the m-th Monte Carlo sample
      mZ[m, ] = parSapply(cl = cluster, 1:iB, function(b, lFoo1){
        set.seed(b)
        return(get_testStat.fast(lFoo1$Y, lFoo1$X, dNu, dTeff))
      }, lFoo1 = lFoo1)
    }

    # Create a vector of Qs across Monte Carlo samples
    vQ = sapply(1:iM, function(m){mean(mZ[m, ]<=dC)})

    # Store rejection frequencies for both thresholds
    mRes1[n, j] = mean(vQ < dCritic1)
    mRes2[n, j] = mean(vQ < dCritic2)
  }
}



stopCluster(cluster)
rm(list=ls())
gc()



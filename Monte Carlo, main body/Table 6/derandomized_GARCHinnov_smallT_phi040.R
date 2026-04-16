library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <- 5 # Numbr of cores
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

# Load the data. Change according to your working directory
sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Size/"
sType = "Size"
load(paste(sPathData, sType, "LatentFactor_PhiNu040_garch.Rdata", sep = ""))

dAlpha = 0.05 #Level of the one shot tests
dNu <- 5 # Power for psi
dQ = 1/4 # Power for f(B)

# Pre-allocation matrices for rejection frequencies using different thresholds
mRes1 = matrix(0.0, ncol = length(vT), nrow = length(vN))
mRes2 = matrix(0.0, ncol = length(vT), nrow = length(vN))

for(n  in 1:length(vN)){
  iN = vN[n]
  iB = (floor(log(iN)^2)) # Number of trials per Monte Carlo samples
  
  mZ = matrix(0.0, nrow = iM, ncol =iB)
  
  # Normalizing sequences for EVT asymptotics and critical value 
  lAB = get_AB(iN)
  dC = get_criticalValues_ultimate(lAB$A, lAB$B, dAlpha)
  
  lFoo = lData[[n]]
  
  # LIL and non-LIL based thresholds
  dCritic1 = (1.0-dAlpha) - sqrt(dAlpha*(1.0-dAlpha))*sqrt(2.0*log(log(iB))/iB)
  dCritic2 = (1.0-dAlpha) - iB^(-dQ)
  
  for(j in 1:length(vT)){
    iT = vT[j]
    dTeff <- iT^(1/dNu)
    
    for(m in 1:iM){
      set.seed(m)
      lFoo1 = list(Y = lFoo[[m]]$Y[, 1:iT], X = lFoo[[m]]$X[, 1:iT])
      clusterExport(cl=cluster, c("dNu", "dTeff"))
      
      # Calculate test statistics across the B trials
      mZ[m, ] = parSapply(cl = cluster, 1:iB, function(b, lFoo1){
        set.seed(b)
        return(get_testStat.fast(lFoo1$Y, lFoo1$X, dNu, dTeff))
      }, lFoo1 = lFoo1)
    }
    
    # Calculate Q across Monte Carlo samples
    vQ = sapply(1:iM, function(m){mean(mZ[m, ]<=dC)})
    
    # Store rejection frequencies
    mRes1[n, j] = mean(vQ < dCritic1)
    mRes2[n, j] = mean(vQ < dCritic2)
    print(mRes1)
    print(mRes2)
  }
}



stopCluster(cluster)
rm(list=ls())
gc()



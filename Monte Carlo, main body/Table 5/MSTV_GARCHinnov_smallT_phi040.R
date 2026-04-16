library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <- 10
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

dAlpha = 0.05 # Significance level
iM = 1000 # Number of MC samples

dNu = 5 # power for psi

vN = c(100, 200, 500) # Cross-sectional sizes
vT = c(100, 200, 300, 500) # Temporal sizes

# Rejection frequencies under the null or the alternative
sType = "size"
# sType = "power"

if(sType == "power"){
  sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Power/"
} else {
  sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Size/"
}
load(paste(sPathData, sType, "LatentFactor_PhiNu040_garch.Rdata", sep = ""))

# Pre-allocation matrix for rejection frequencies
mRej = matrix(0.0, ncol = length(vT), nrow = length(vN))


for(n  in 1:length(vN)){
  iN = vN[n]
  lFoo = lData[[n]]
  lAB = get_AB(iN) # Normalizing sequences for EVT asymptotics
  dC = get_criticalValues_ultimate(lAB$A, lAB$B, dAlpha) # critical values

  for(j in 1:length(vT)){
    iT = vT[j]
    dTeff = iT^(1/dNu)
    clusterExport(cl =cluster, c("dNu", "dTeff", "iT"))

    # Calculate the statistic
    vZ = parSapply(cl = cluster, seq(iM), function(m, lFoo){
      set.seed(m)
      return(get_testStat.fast(lFoo[[m]]$Y[, 1:iT], lFoo[[m]]$X[, 1:iT], dNu, dTeff))
    }, lFoo = lFoo)

    # Store rejectionf requencies
    mRej[n,j] = mean(vZ > dC)
  }
}

mRej

stopCluster(cluster)
rm(list=ls())
gc()



library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <- 10 # Number of cores
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
vN = c(100, 200, 500) #Cross-sectional sizes
vT = c(100, 200, 300, 400, 500) # Temporal sizes


# Repository, change according to your PC
sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Size/"
sType = "Size"


dC = qcauchy(0.95) # Critical value
mRej = matrix(0.0, ncol = length(vT), nrow = length(vN)) # Pre-allocation matrix
load(paste(sPathData, sType, "LatentFactor_PhiNu040_T.Rdata", sep = ""))


for(n  in 1:length(vN)){
  iN = vN[n]
  lFoo = lData[[n]]

  for(j in 1:length(vT)){
    iT = vT[j]

    # Calculate test statistic
    vTstat = parSapply(cl=cluster, seq(iM), function(m, lFoo, iT){
      cval <- SCTcv(t(lFoo[[m]]$Y[, 1:iT]), t(lFoo[[m]]$X[, 1:iT]), iL =0)
      return(cval)
    }, lFoo=lFoo, iT = iT)
    # Store rejection frequencies
    mRej[n,j] = mean(vTstat > dC)
  }
}

stopCluster(cluster)
rm(list=ls())
gc()

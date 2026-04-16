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
iM = 1000 # Number of MC samples

vN = c(100, 200, 500) # Cross-sectional sizes
vT = c(100, 200, 300, 500) # Temporal sizes

# Path for the simulate data, change according to your settings
sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Size/"
sType = "Size"
load(paste(sPathData, sType, "LatentFactor_PhiNu040_T.Rdata", sep = ""))

dC = get_criticalValuesFLLM(0.05) # Critical values
mRej = matrix(0.0, ncol = length(vT), nrow = length(vN)) # pre-allocation matrix

for(n  in 1:length(vN)){
  iN = vN[n]
  lFoo = lData[[n]]

  for(j in 1:length(vT)){
    iT = vT[j]

    # Calculate test statistics
    vZ = parSapply(cl=cluster, 1:iM, function(m, lFoo, iT){
          get_testStatFLLM(lFoo[[m]]$Y[, 1:iT], lFoo[[m]]$X[, 1:iT])
        }, lFoo = lFoo, iT = iT)

    # Store rejection frequencies
    mRej[n,j] = mean(vZ  - 2.0*log(iN) + log(log(iN))> dC)
  }
}

mRej # print matrix of rejection frequencies

stopCluster(cluster)
rm(list=ls())
gc()


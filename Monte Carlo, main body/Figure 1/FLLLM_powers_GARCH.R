library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <-8 # Number of cores
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

iK = 3 # Number of factors
iM = 1000 # Number of MC samples

# Parameters for Factors DGP 
vIntF = c(0.53, 0.19, 0.19)
mPhi = diag(c(-0.1, 0.2, -0.2))


# Critical value
dAlpha = 0.05
dC = get_criticalValuesFLLM(dAlpha)


# Possible values of \vartheta
boo = seq(1,5, 0.5)
# Pre-allocation vector for rejection frequencies
vRej = numeric(length(boo))
iN = 500 # cross-sectional size
iT = 100 # temporal size
clusterExport(cl =cluster, "iT")


for(k in 1:length(boo)){
  # Simulate data
  lInput = list(iK = iK, iN = iN, iT = iT, Phi = mPhi, Type="power",
                PhiNu = 0.40, vIntF = vIntF, dPerc = 0.95, boo = boo[k])
  lData = get_data_latentFactor_GARCH.bis_blocks2(lInput)
  
  print(paste("Data for boo =", boo[k], " generated"))
  clusterExport(cl =cluster, c("lData"))
  
  # Calculate test statistics across MC samples
  vZ = parSapply(cl = cluster, seq(iM), function(m){
    set.seed(m)
    return(get_testStatFLLM(lData[[m]]$Y[, 1:iT], lData[[m]]$X[, 1:iT]))
  })
  
  # Store rejection frequencies
  vRej[k] = mean(vZ  - 2.0*log(iN) + log(log(iN)) > dC)
}


stopCluster(cluster)
rm(list=ls())
gc()

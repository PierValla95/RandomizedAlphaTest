
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

sType = "Power"
dC = qnorm(.95) # Critical value

# Possible values of \varheta
boo = seq(1,5, 0.5)

# Pre-allocation vector for rejection frequencies
vRej = numeric(length(boo))

iN = 500 # Cross-sectional sizes
iT = 100 # Temporal sizes
clusterExport(cl =cluster, "iT")

# Loop across possible avlues of \vartheta
for(k in 1:length(boo)){
  
  # Simulate data
  lInput = list(iK = iK, iN = iN, iT = iT, Phi = mPhi, Type="power",
                PhiNu = 0.40, vIntF = vIntF, dPerc = 0.95, boo = boo[k])
  lData = get_data_latentFactor_GARCH.bis_blocks2(lInput)
  clusterExport(cl =cluster, c("lData"))
  
  # Calculate test statistics across MC samples
  vZ = parSapply(cl = cluster, seq(iM), function(m){
    set.seed(m)
    return(get_testStatPY(lData[[m]]$Y[, 1:iT], lData[[m]]$X[, 1:iT]))
  })
   
  # Store rejection frequencies
  vRej[k] = mean(vZ > dC)
}




stopCluster(cluster)
rm(list=ls())
gc()

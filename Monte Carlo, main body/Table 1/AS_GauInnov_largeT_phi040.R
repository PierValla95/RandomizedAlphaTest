
library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <-5
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

vN = c(100, 200, 500) # Cross-sectional sized
vT = c(1000, 2000) # Temporal sizes
iT.max = max(vT)

dC = qcauchy(.95) # Critical value

# Pre-allocation matrix for rejection frequencies
mRej = matrix(0.0, ncol = length(vT), nrow = length(vN))
for(n in seq_along(vN)){
  iN = vN[n]

  # Simulate data. Set Type = "size" for the null, Type = "power" for the alternative
  lInput = list(iK = iK, iN = iN, iT = iT.max, Phi = mPhi, Type="size",
                PhiNu = 0.40, vIntF = vIntF, dPerc = 0.95)
  lData = get_data_latentFactor.bis(lInput)
  clusterExport(cl =cluster,"lData")
  
  for(j in seq_along(vT)){
    iT = vT[j]
    clusterExport(cl =cluster, "iT")
    
    # Calculate the test statistics
    vTstat = parSapply(cl=cluster, seq(iM), function(m){
      cval <- SCTcv(t(lData[[m]]$Y[, 1:iT]), t(lData[[m]]$X[, 1:iT]), iL =0)
      return(cval)
    })

    # Store rejection frequencies
    mRej[n,j] = mean(vTstat > dC)
  }
}

stopCluster(cluster)
rm(list=ls())
gc()


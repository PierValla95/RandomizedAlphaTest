library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <-210 # Number of cores
set.seed(150)
cluster = makeCluster(iH)
clusterEvalQ(cluster, {
  library(RandAlphaTest)
})

iK = 3 # Number of factors
iM = 1000 # Number of MC samples

# Parameters for Factors DGP
vIntF = c(0.53, 0.19, 0.19)
mPhi = diag(c(-0.1, 0.2, -0.2))

vN = c(100, 200, 500) # Cross sectional sizes
vT = c(1000, 2000) # Temporal sizes
iT.max = max(vT)

dNu = 5 # Exponent of psi
dAlpha = 0.05 # Level of the test

# Pre-allocaton matrix
mRej = matrix(0.0, ncol = length(vT), nrow = length(vN))

for(n in seq_along(vN)){
  iN = vN[n]

  # Simulate data. Choose "size" for under the null and "power" for the alternative
  lInput = list(iK = iK, iN = iN, iT = iT.max, Phi = mPhi, Type="size",
                PhiNu = 0.40, vIntF = vIntF, dPerc = 0.95)
  lData = get_data_latentFactor.bis(lInput)

  # Normalizing sequences and critical value
  lAB = get_AB(iN)
  dC = get_criticalValues_ultimate(lAB$A, lAB$B, dAlpha)
  clusterExport(cl =cluster, c("dNu", "lData"))

  for(j in seq_along(vT)){
    iT = vT[j]
    dTeff = iT^(1/dNu)
    clusterExport(cl =cluster, c("dTeff", "iT"))

    # Calculate test statistics
    vZ = parSapply(cl = cluster, seq(iM), function(m){
      set.seed(m)
      return(get_testStat.fast(lData[[m]]$Y[, 1:iT], lData[[m]]$X[, 1:iT], dNu, dTeff))
    })

    # Store rejection frequencies
    mRej[n,j] = mean(vZ > dC)
  }
}

stopCluster(cluster)
rm(list=ls())
gc()

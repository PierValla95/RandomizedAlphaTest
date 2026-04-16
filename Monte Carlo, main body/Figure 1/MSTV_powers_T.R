library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <-8
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

dNu = 5 # Power of psi
dAlpha = 0.05 # Level of the test

boo = seq(1, 5, 0.50) # Possible values of \vartheta 

vRej = numeric(length(boo)) # Pre-allocation vector for rejection frequencies
iN = 500 # Cross-sectional size

# Normalizing sequences for EVT asymptotics and critical value
lAB = get_AB(iN)
dC = get_criticalValues_ultimate(lAB$A, lAB$B, dAlpha)

iT = 100 # Temporal size
dTeff = iT^(1/dNu)
clusterExport(cl =cluster, c("dTeff", "iT"))

# Loop across possible values of \vartheta
for(k in 1:length(boo)){
  
  # Simulate data
  lInput = list(iK = iK, iN = iN, iT = iT, Phi = mPhi, Type="power",
                PhiNu = 0.40, vIntF = vIntF, dPerc = 0.95, boo = boo[k])
  lData = get_data_latentFactor_T.bis_blocks2(lInput)
  
  clusterExport(cl =cluster, c("dNu", "lData"))
  
  # Calculate test statistics across MC samples
  vZ = parSapply(cl = cluster, seq(iM), function(m){
    set.seed(m)
    return(get_testStat.fast(lData[[m]]$Y[, 1:iT], lData[[m]]$X[, 1:iT], dNu, dTeff))
  })
  vRej[k] = mean(vZ > dC)
}


stopCluster(cluster)
rm(list=ls())
gc()

library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <-5 # Number of cores
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

vN = c(100, 200, 500) # Cross-sectional sizes
vT = c(1000, 2000) # Temporal sizes
iT.max = max(vT)


dQ = 1/4 # Power for f(B)
dNu = 5 # Power for psi
dAlpha = 0.05 # Level of the one-shot tests


# Pre-allocation matrices for rejection frequencies
mRes1 = matrix(0.0, ncol = length(vT), nrow = length(vN))
mRes2 = matrix(0.0, ncol = length(vT), nrow = length(vN))

for(n in seq_along(vN)){
  iN = vN[n]
  
  # Simulate the data. Set Type = "size" for under the null and Type = "power" for the alternative
  lInput = list(iK = iK, iN = iN, iT = iT.max, Phi = mPhi, Type="size",
                PhiNu = 0.40, vIntF = vIntF, dPerc = 0.95)
  lData = get_data_latentFactor_withGarch.bis(lInput)
  
  # Normalizing sequences for EVT asymptotics and critical value
  lAB = get_AB(iN)
  dC = get_criticalValues_ultimate(lAB$A, lAB$B, dAlpha)
  
  # Number of trials for each MC sample
  iB = (floor(log(iN)^2))
  
  mZ = matrix(0.0, nrow = iM, ncol =iB)
  
  # LIL and non-LIL based thresholds
  dCritic1 = (1.0-dAlpha) - sqrt(dAlpha*(1.0-dAlpha))*sqrt(2.0*log(log(iB))/iB)
  dCritic2 = (1.0-dAlpha) - iB^(-dQ)
  
  
  for(j in seq_along(vT)){
    iT = vT[j]
    dTeff = iT^(1/dNu)
    vIdx = seq(iT)
    for(m in seq(iM)){
      set.seed(m)
      
      # Calculate the statstics Psi only once, not for all B trials
      mY = lData[[m]]$Y[, vIdx]
      mX = cbind(rep(1, iT), t(lData[[m]]$X[, vIdx]))
      mFoo.x = solve(crossprod(mX))
      mOLS = mFoo.x %*% crossprod(mX, t(mY))
      vAlpha = mOLS[1, ]
      mRes = t(mY) - mX %*% mOLS
      vPsi = (dTeff * (abs(vAlpha)/sqrt(mean(mRes^2))))^(dNu/2)
      
      # Calculate the block maximum for each trial
      mZ[m, ] = parSapply(cl = cluster, seq(iB), function(b, vPsi, iN){
        set.seed(b)
        return(max(vPsi + rnorm(iN)))
      }, vPsi = vPsi, iN = iN)}
    
    # Calculate Q across the MC samples
    vQ = sapply(1:iM, function(m){mean(mZ[m, ]<=dC)})
    
    # Store rejection frequencies
    mRes1[n, j] = mean(vQ < dCritic1)
    mRes2[n, j] = mean(vQ < dCritic2)
  }
}

stopCluster(cluster)
rm(list=ls())
gc()


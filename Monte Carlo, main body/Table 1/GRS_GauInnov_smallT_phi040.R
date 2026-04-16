rm(list=ls())
gc()
library(RandAlphaTest)

dAlpha = 0.05 # Significance level
iM = 1000 # Number of MC samples

vN = c(100, 200, 500) # Cross-sectional sizes
vT = c(100, 200, 300, 500) # Temporal sizes

sType = "Size" # Size or power
# Repository, change according to your PC
sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Size/"
# Load the data
load(paste(sPathData, sType, "LatentFactor_PhiNu040.Rdata", sep = ""))

# Pre-allocation matrix for rejection frequencies
mRej = matrix(0.0, ncol = length(vT), nrow = length(vN))

for(n  in 1:length(vN)){
  iN = vN[n]
  lFoo = lData[[n]]

  for(j in 1:length(vT)){
    iT = vT[j]
    # GRS can only be calculated when T > N
    if(iT > iN){
      dC = qf(0.95, df1 = iT - iN - 3, df2 = iN) # Critical values

      # Calculate test statistic
      vZ = sapply(1:iM, function(m){
        get_testStatGRS(lFoo[[m]]$Y[, 1:iT], lFoo[[m]]$X[, 1:iT], iT, iN)
      })

      # Store rejection frequencies
      mRej[n,j] = mean(vZ > dC)
    } else {
      mRej[n,j] = 0.0
    }
  }
}



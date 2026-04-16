library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <- 15
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


vN = c(100, 200, 500)
vT = c(100, 200, 300, 500)

sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Power/"

dAlpha = 0.05

dNu <- 5
dQ = 1/4

sType = "power"
mRes1 = matrix(0.0, ncol = length(vT), nrow = length(vN))
mRes2 = matrix(0.0, ncol = length(vT), nrow = length(vN))
load(paste(sPathData, sType, "LatentFactor_PhiNu040.Rdata", sep = ""))
for(n  in 1:length(vN)){
  iN = vN[n]
  iB = (floor(log(iN)^2))

  mZ = matrix(0.0, nrow = iM, ncol =iB)

  lAB = get_AB(iN)
  dC = get_criticalValues_ultimate(lAB$A, lAB$B, dAlpha)
  lFoo = lData[[n]]

  dCritic1 = (1.0-dAlpha) - sqrt(dAlpha*(1.0-dAlpha))*sqrt(2.0*log(log(iB))/iB)
  dCritic2 = (1.0-dAlpha) - iB^(-dQ)

  for(j in 1:length(vT)){
    iT = vT[j]
    dTeff <- iT^(1/dNu)

    for(m in 1:iM){
      set.seed(m)
      lFoo1 = list(Y = lFoo[[m]]$Y[, 1:iT], X = lFoo[[m]]$X[, 1:iT])
      clusterExport(cl=cluster, c("dNu", "dTeff"))
      mZ[m, ] = parSapply(cl = cluster, 1:iB, function(b, lFoo1){
        set.seed(b)
        return(get_testStat.fast(lFoo1$Y, lFoo1$X, dNu, dTeff))
      }, lFoo1 = lFoo1)
    }

    vQ = sapply(1:iM, function(m){mean(mZ[m, ]<=dC)})
    mRes1[n, j] = mean(vQ < dCritic1)
    mRes2[n, j] = mean(vQ < dCritic2)
    print(c(iT, iN, mRes1[n, j], mRes2[n, j]))
  }
}



stopCluster(cluster)
rm(list=ls())
gc()



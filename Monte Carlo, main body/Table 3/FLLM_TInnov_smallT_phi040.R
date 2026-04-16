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

vN = c(100, 200, 500)
vT = c(100, 200, 300, 500)

sPathData = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Size/"
sType = "Size"

dC = get_criticalValuesFLLM(0.05)
mRej = matrix(0.0, ncol = length(vT), nrow = length(vN))
load(paste(sPathData, sType, "LatentFactor_PhiNu040_T.Rdata", sep = ""))

for(n  in 1:length(vN)){
  iN = vN[n]
  lFoo = lData[[n]]

  for(j in 1:length(vT)){
    iT = vT[j]

    vZ = parSapply(cl=cluster, 1:iM, function(m, lFoo, iT){
      get_testStatFLLM(lFoo[[m]]$Y[, 1:iT], lFoo[[m]]$X[, 1:iT])
    }, lFoo = lFoo, iT = iT)


    mRej[n,j] = mean(vZ  - 2.0*log(iN) + log(log(iN))> dC)
    print(c(iN, iT, mRej[n,j]))
  }
}


stopCluster(cluster)
rm(list=ls())
gc()


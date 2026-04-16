## THIS SCRIPT SIMULATE DATA WITH STUDENT'S T INNOVATIONS FROM MASSACCI ET AL. (2026)

library(RandAlphaTest)
library(parallel)

# For parallel computing
iH <-10 # Number of cores
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
vT = c(100, 200, 300, 500) # Temporal sizes
iT.max = max(vT)

lData = list() # Pre-allocation list

for(n in seq_along(vN)){
  iN = vN[n]

  # Set Type = "size" for data under the null, Type = "power" for data under the alternative
  # Set "PhiNu = 0.00" for no persistence in the omitted factor, "PhiNu = 0.40" for persistence as in the main body
  lInput = list(iK = iK, iN = iN, iT = iT.max, Phi = mPhi, Type="size",
                PhiNu = 0.40, vIntF = vIntF, dPerc = 0.95)
  lFoo = get_data_latentFactor_T.bis(lInput)
  lData[[n]] = lFoo
}

# Save data in your preferred working directory.
# Use "Size" at the beginning of the name for data under the null
# Use "Power" at the beginning of the name for data under the null
# USe "PhiNu040" for data with persistence in the omitted factor (Monte Carlo of the main body)
# USe "PhiNu00" for data without  persistence in the omitted factor (Monte Carlo of the supplement)

save(lData, file = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Simulations/sim_data/Size/SizeLatentFactor_PhiNu040_T.Rdata")

stopCluster(cluster)
rm(list=ls())
gc()


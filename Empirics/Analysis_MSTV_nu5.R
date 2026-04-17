library("RandAlphaTest")
library("parallel")


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

# Replace with appropriate path to data
sPath = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Empirics/Data/"

# Load returns on cross-section of S&P500 constituents
filename="data_S&P500_fullSample_min5Yrs.Rdata"
load(file = paste(sPath, filename, sep =""))
mRet = lData$ret
vDates = lData$dates

# Load FF Factors
mFF = read.csv(paste(sPath, "FF_Data_new.csv", sep = ""), sep = ";")

# Align dates
iIdx.low = which(vDates == "1985-01-01")
mRet = mRet[(iIdx.low):nrow(mRet), ]
mFF = mFF[(iIdx.low):nrow(mFF), ]

iR = 60 # Length of the window
iK = nrow(mRet) - iR +1
dNu = 5 # Value of nu in the test statistics
dAlpha = 0.05 #Significance level

# Pre allocation
model_names <- c("capm", "FF2", "FF3", "FF4", "FF5", "FF6")
lRes <- lQ <- setNames(replicate(length(model_names), numeric(iK), simplify = FALSE), model_names)
vCritic2 <- numeric(iK)

dQ = 1/4 # Power for the non-lilthreshold in the derandomized procedure

# Possible models
model_factors <- list(capm = 1, FF2 = c(1, 7), FF3 = 1:3, FF4 = c(1:3, 7),
                      FF5 = 1:5, FF6 = c(1:5, 7))


for(k in 1:iK){
  print(round(100*k/iK, 3))
  set.seed(k) #set seed for replicability

  # Consturct the window
  iL = k
  iU = iL + iR - 1
  mRet1 = mRet[iL:iU, ]
  mFactors = mFF[iL:iU, -1] # All factors + RF


  # Remove returns above 200%. This happens for one/two stocks in less than 80 windows.
  mRet1[!is.na(mRet1) & abs(mRet1) > 200] <- NA

  # Remove missing values
  vIdx = which(colSums(is.na(mRet1)) > 0)
  if(sum(vIdx)>0){
    mRet1 = mRet1[, -vIdx]
  }

  iN = ncol(mRet1) # Number of assets
  vRF <- mFactors[, 6] # Risk-free rate
  mY <- mRet1 - matrix(rep(vRF, iN), nrow = iR, ncol = iN) # Excess returns


  iB = (floor(log(iN)^2)) # Number of repetitions for the derandomized procedure
  vCritic2[k] =  (1.0-dAlpha) - iB^(-dQ) #Critical values for de-randomized procedure

  lAB = get_AB(iN) # Norming sequences for EVT asymptotics
  dC = get_criticalValues_ultimate(lAB$A, lAB$B, dAlpha) # Critical values based on Gumbel quantiles

  # Loop the procedure across models
  for (mod in model_names) {
    X <- t(mFactors[, model_factors[[mod]]])
    Y <- t(mY)
    lInput <- list(Y = Y, X = X, nu = dNu)

    # Calculate the test statistics across B trials
    vZ <- parLapply(cl = cluster, 1:iB, function(b, lInput) {
      set.seed(b)
      get_testStat_unbalanced(lInput$Y, lInput$X, lInput$nu)
    }, lInput = lInput)

    # Rejection results
    lQ[[mod]][k] <- mean(unlist(vZ) <= dC)
    lRes[[mod]][k] <- lQ[[mod]][k] < vCritic2[k]
  }
}

stopCluster(cluster)


# Combine results
mRes <- do.call(cbind, lRes)
round(colMeans(mRes), 3) # Values in the upper panel of Table 5.1

# Indexes for period of market downturn
vRecL = c(91, 123, 216, 361)
vRecU = c(108, 154, 234, 377)

# Crisis-specific rejections
for(i in seq_along(vRecL)){
  print(round(colMeans(mRes[vRecL[i]:vRecU[i], ]),3))
}


# Rejections over all crisis and non-crisis periods
vIdx = c(91:108, 123:154, 216:234, 361:377)
round(colMeans(mRes[vIdx, ]), 3)
round(colMeans(mRes[-vIdx, ]), 3)











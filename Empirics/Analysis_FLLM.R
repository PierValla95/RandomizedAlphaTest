library("RandAlphaTest")


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
iK = nrow(mRet) - iR +1 # Number of windows
dAlpha = 0.05 #Significance level

# Pre allocsation
model_names <- c("capm", "FF2", "FF3", "FF4", "FF5", "FF6")
lRes <- setNames(replicate(length(model_names), numeric(iK), simplify = FALSE), model_names)

model_factors <- list(capm = 1, FF2 = c(1, 7), FF3 = 1:3, FF4 = c(1:3, 7),
                      FF5 = 1:5, FF6 = c(1:5, 7))

# Critical values
dC = get_criticalValuesFLLM(dAlpha)


for(k in 1:iK){
  print(round(100*k/iK, 3)) # Counter

  # Beginning and end of the window
  iL = k
  iU = iL + iR - 1
  mRet1 = mRet[iL:iU, ] # Returns for the window
  mFactors = mFF[iL:iU, -1] # All factors + Risk free rate

  mRet1[!is.na(mRet1) & abs(mRet1) > 200] <- NA # Remove NA and returns above 200% (this happens in few windows for one or two assets out of more than 500)
  vIdx = which(colSums(is.na(mRet1)) > 0)
  if(sum(vIdx)>0){
    mRet1 = mRet1[, -vIdx] # No missing data
  }

  iN = ncol(mRet1) # Number of assets
  vRF <- mFactors[, 6] # Risk free rate
  mY <- mRet1 - matrix(rep(vRF, iN), nrow = iR, ncol = iN) # Calculate excess returns

  for (mod in model_names) {
    dZ = get_testStatFLLM(t(mY),  t(mFactors[, model_factors[[mod]]]))
    lRes[[mod]][k] <- ((dZ - 2.0*log(iN) + log(log(iN))) > dC)
  }
}

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






library("RandAlphaTest")

# Load test assets returns
sPath = "C:/Users/pierl/Dropbox/Massacci_Sarno_Trapani2/Empirics/Data/"
filename="data_S&P500_fullSample_min5Yrs.Rdata"
load(file = paste(sPath, filename, sep =""))
mRet = lData$ret
vDates = lData$dates

# Load Fama French factors
mFF = read.csv(paste(sPath, "FF_Data_new.csv", sep = ""), sep = ";")

# Prepare the data
iIdx.low = which(vDates == "1985-01-01")
mRet = mRet[(iIdx.low):nrow(mRet), ]
mFF = mFF[(iIdx.low):nrow(mFF), ]

iR = 60 # Length of the window
iK = nrow(mRet) - iR +1 # Number of windows

# Possible models
model_factors <- list(capm = 1, FF2 = c(1, 7), FF3 = 1:3, FF4 = c(1:3, 7),
                      FF5 = 1:5, FF6 = c(1:5, 7))

# Pre-allocation lists
lAlpha = list()
lVar.res = list()

mod = "FF2" # Focus on CAPM+Momentum

 # Loop across windows
for(k in 1:iK){
  set.seed(k) #set seed for replicability

  # Construct the window
  iL = k
  iU = iL + iR - 1
  mRet1 = mRet[iL:iU, ]
  mFactors = mFF[iL:iU, -1] # All factors + RF

   # Remove missing entries and returns above 200% in one month
  mRet1[!is.na(mRet1) & abs(mRet1) > 200] <- NA
  vIdx = which(colSums(is.na(mRet1)) > 0)
  if(sum(vIdx)>0){
    mRet1 = mRet1[, -vIdx] # No missing data
  }

  iN = ncol(mRet1) # Number of assets
  vRF <- mFactors[, 6] # Risk-free rate

  # excess returns and factors
  mY <- mRet1 - matrix(rep(vRF, iN), nrow = iR, ncol = iN)
  mF = t(mFactors[, model_factors[[mod]]])

  # OLS estimation
  mX = cbind(rep(1.0, iR), t(mF))
  mFoo.x =  (solve(crossprod(mX)))
  mOLS = mFoo.x%*%crossprod(mX, mY)

  # Store vector of alphas
  lAlpha[[k]] = mOLS[1,]

  # Store residuals' variance
  mRes = mY - mX%*%mOLS
  lVar.res[[k]] = apply(mRes, 2, var)
}



# Dates for the plot
vDates.plot =  seq(from = as.Date("1990-01-01"),
                   to = as.Date("2024-12-01"),
                   by = "month")

# Indexes of recession dates
vRecL = c(91, 123, 216, 361)
vRecU = c(108, 154, 234, 377)

# Ticks for the x-axis and labels
x_ticks = as.Date(c("1992-01-01", "1995-01-01", "1997-01-01", "2000-01-01",
                    "2002-01-01", "2005-01-01", "2007-01-01", "2010-01-01",
                    "2012-01-01", "2015-01-01", "2017-01-01", "2020-01-01",
                    "2022-01-01", "2025-01-01"))
x_labels = format(x_ticks, "%Y")

# Indices where FLLM rejects very few times while ours does
iIdx.low = which(vDates.plot == "1999-01-01")
iIdx.Up = which(vDates.plot == "2009-06-01")

# Five largest absolute pricing errors over the K widnows
mAlpha.fiveLargest =do.call(rbind, lapply(1:iK, function(k){
  tail(sort(abs(lAlpha[[k]])), 5)
}))

# Residuals' variance of the most misspriced assets
mVar.fiveLargest =do.call(rbind, lapply(1:iK, function(k){
  foo = tail(order(abs(lAlpha[[k]])), 5)
  return(lVar.res[[k]][foo])
}))


# Plot the time series for Figure 5.1
k = 0
for(i in 5:3){
  k = k+1
  windows(5, 5)
  plot(vDates.plot, mAlpha.fiveLargest[, i],  xaxt = "n", type = "n", lwd = 2, las = 1,
       xlab = "", ylab = "")
  axis(1, at = x_ticks, labels = x_labels)
  grid()
  vX1 = vDates.plot[iIdx.low:iIdx.Up]
  polygon(x = c(vX1, rev(vX1)), y =  c(rep(0.0, length(vX1)),
                                       rev(rep(1.1*max(mAlpha.fiveLargest[, 5]), length(vX1)))),
          col =  adjustcolor("red", alpha.f = 0.3), border = NA)
  lines(vDates.plot,  mAlpha.fiveLargest[, i], lwd = 2, col = 1)


  windows(5, 5)
  plot(vDates.plot, sqrt(mVar.fiveLargest[, i]), xaxt = "n", type = "n", lwd = 2, las = 1,,
       xlab = "", ylab = "")
  axis(1, at = x_ticks, labels = x_labels)
  grid()
  vX1 = vDates.plot[iIdx.low:iIdx.Up]
  polygon(x = c(vX1, rev(vX1)), y =  c(rep(0.0, length(vX1)),
                                       rev(rep(1.1*max(mVar.fiveLargest[, 5]), length(vX1)))),
          col =  adjustcolor("red", alpha.f = 0.3), border = NA)
  lines(vDates.plot,  sqrt(mVar.fiveLargest[, i]), lwd = 2, col = 1)
}






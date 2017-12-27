# Calculate VaR and ES at 5% level

# Size of train sample
N <- round(n*0.25)
testRet <- xts(tail(portRet, -N), date[(N+1):n])
testDate <- date[(N+1):n]

# Empirical
empVaR.curve <- empES.curve <- xts(rep(0, n-N), testDate)
for (i in 1:(n-N)) {
	trainRet <- portRet[i:(i+N-1)]
	empVaR.curve[i] <- -VaR(trainRet)#, method = 'historical')
	empES.curve[i] <- -ES(trainRet, method='historical')
}
# Kupiec test
empVaR.test <- kupiec.test(testRet, empVaR.curve)

# Gaussian copula
nVaR.curve <- nES.curve <- xts(rep(0, n-N), testDate)
for (i in 1:(n-N)) {
	trainRet <- as.matrix(log_returns[i:(i+N-1)])
	trainData <- apply(trainRet, 2, function(p) quantile(p, nSimCop))
	trainCopRet <- as.numeric(trainData %*% weights)
	nVaR.curve[i] <- -VaR(trainCopRet)#, method = 'historical')
	nES.curve[i] <- -ES(trainCopRet, method='historical')
	cat('=')
}
# Kupiec test
nVaR.test <- kupiec.test(testRet, nVaR.curve)

# Gaussian copula
tVaR.curve <- tES.curve <- xts(rep(0, n-N), testDate)
for (i in 1:(n-N)) {
	trainRet <- as.matrix(log_returns[i:(i+N-1)])
	trainData <- apply(trainRet, 2, function(p) quantile(p, tSimCop))
	trainCopRet <- as.numeric(trainData %*% weights)
	tVaR.curve[i] <- -VaR(trainCopRet)#, method = 'historical')
	tES.curve[i] <- -ES(trainCopRet, method='historical')
	cat('=')
}
# Kupiec test
tVaR.test <- kupiec.test(testRet, tVaR.curve)

# Vine copula
vcVaR.curve <- vcES.curve <- xts(rep(0, n-N), testDate)
for (i in 1:(n-N)) {
	trainRet <- as.matrix(log_returns[i:(i+N-1)])
	trainData <- apply(trainRet, 2, function(p) quantile(p, simCop))
	trainCopRet <- as.numeric(trainData %*% weights)
	vcVaR.curve[i] <- -VaR(trainCopRet)#, method = 'historical')
	vcES.curve[i] <- -ES(trainCopRet, method='historical')
	cat('=')
}
# Kupiec test
vcVaR.test <- kupiec.test(testRet, vcVaR.curve)

plotVaR <- merge.xts(testRet, -empVaR.curve, -nVaR.curve, -tVaR.curve, -vcVaR.curve)
plotES <- merge.xts(testRet, -empES.curve, -nES.curve, -tES.curve, -vcES.curve)
colnames(plotVaR) <- colnames(plotES) <- c('Profit & Loss', 'Empirical', 'Gaussian copula', "Student's t copula", 'Vine copula')


# Summary plot
pdf(file='Export/VaRcurve.pdf', width=10, height=6)
# setEPS(); postscript("Export/VaRcurve.eps")
chart.TimeSeries(plotVaR, xlab='Date', ylab='Profit & Loss', main="", date.format = "%d.%m.%Y", 
                 lty=c(1,2,1,1,2), colorset=c('black', 'black', 'blue', 'red', 'violet'), 
                 legend.loc='topleft',  pch=NULL, lwd=1)
# "Profit & Loss and 5% VaR curve"
dev.off()

pdf(file='Export/EScurve.pdf', width=10, height=6)
# setEPS(); postscript("Export/EScurve.eps")
chart.TimeSeries(plotES, main="", xlab='Date', ylab='Profit & Loss', date.format = "%d.%m.%Y", 
                 lty=c(1,2,1,1,2), colorset=c('black', 'black', 'blue', 'red', 'violet'), 
                 legend.loc='topleft',  pch=NULL, lwd=1)
# "Profit & Loss and 5% ES curve"
dev.off()
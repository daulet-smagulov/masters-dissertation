# rm(list=ls(all=TRUE))
# setwd("/Users/Daulet/Documents/University/Masters/Science/SP/Main")

# Goodness of fit tests
if(!require(goftest)) {install.packages("goftest"); library(goftest)}

source("Initialization.r")

x <- seq(min(log_returns), max(log_returns), length.out=100)

################# HYPERBOLIC DISTRIBUTION #################

# # Simulate theorerical data
hSim <- sapply(size, function(p) rhyperb(n, hTheta[,p]))
colnames(hSim) = paste("hSim for",assets)
hSimPDF <- sapply(size, function(p) dhyperb(hSim[,p], hTheta[,p]))
colnames(hSimPDF) = paste("hSimPDF for",assets)

# # Plot results
# pdf(file='Export/Fitting hyperbolic distribution.pdf', encoding='KOI8-R.enc', width=12, height=7)
setEPS(); postscript("Export/Fitting hyperbolic distribution.eps")
for (i in 1:l) {
	par(mfrow = c(1,2))
	# Show emp histogram with theoretical data line
	hist(log_returns[,i], probability = TRUE, label=TRUE, xlab = paste("Log returns", assets[i]), ylim = c(0,1.1*max(hSimPDF[,i])), main = paste(assets[i], "log returns histogram\nwith theoretical PDF"))
	rug(log_returns[,i])
	lines(x, dhyperb(x, hTheta[,i]), col="red")
	grid()
	# Show Q-Q plot of log returns and simulated data
	qqplot(hSim[,i], as.numeric(log_returns[,i]), main = paste("Q-Q plot of simulated data with\n", assets[i], "log returns"), pch=20, xlab = paste("Log returns", assets[i]), ylab = "Simulated data");
	abline(a=0,b=1, lty=3)
	q05 = qhyperb(0.05, hTheta[,i])
	q95 = qhyperb(0.95, hTheta[,i])
	abline(v=c(q05, 1), lty=2, col='red')
	abline(v=c(q95, 1), lty=2, col='red')
	grid()
}
dev.off()

# # Making two-sample Kolmogorov-Smirnov test
hKSTest = lapply(size, function(p) ks.test(as.numeric(log_returns[,p]), 'phyperb', hTheta[,p]))
for(i in size) {
	hKSTest[[i]]$data.name <- paste("Log returns", assets[i], "and hyperbolic distribution");
	hKSTest[[i]]$method = "One-sample Kolmogorov-Smirnov test of hyperbolic distribution"; }
hKSTest
# hKSTest[[1]]

# Making Anderson-Darling test
hADtest <- lapply(size, function(p) ad.test(log_returns[,p], phyperb, hTheta[,p]))
for(i in size) {
	hADtest[[i]]$data.name <- paste("Log returns", assets[i], "and hyperbolic distr");
	hADtest[[i]]$method = "Anderson-Darling GoF test of hyperbolic distribution"; }
hADtest

# Perform Crämer-von~Mises test
hCvMtest <- lapply(1:l, function(k) cvm.test(log_returns[,k], 'phyperb', hTheta[,k]))
for(i in 1:l) {
	hCvMtest[[i]]$data.name <- paste("Log returns", assets[i])
	hCvMtest[[i]]$method = c("Crämer-von~Mises of GoF", "Null hypothesis: hyperbolic distribution") }
hCvMtest

################### STABLE DISTRIBUTION ###################

# # Simulate theorerical data
stSim = sapply(size, function(i) rstable(n, stAlpha[i], stBeta[i], stGamma[i], stDelta[i])); colnames(stSim) = paste("stSim for",assets)
stSimPDF <- sapply(size, function(i) dstable(stSim[,i], stAlpha[i], stBeta[i], stGamma[i], stDelta[i]))
colnames(stSimPDF) = paste("stSimPDF for",assets)

# # Show emp histogram with theoretical data line
# pdf(file='Export/Fitting stable distribution.pdf', encoding='KOI8-R.enc', width=12, height=7)
setEPS(); postscript("Export/Fitting stable distribution.eps")
for (i in 1:l) {
	par(mfrow = c(1,2))
	# Show emp histogram with theoretical data line
	hist(log_returns[,i], probability = TRUE, label=TRUE, xlab = paste("Log returns", assets[i]), ylim = c(0,1.1*max(stSimPDF[,i])), main = paste(assets[i], "log returns histogram\nwith theoretical PDF"))
	rug(log_returns[,i])
	lines(x, dstable(x, stAlpha[i], stBeta[i], stGamma[i], stDelta[i]), col="red")
	grid()
	# Show Q-Q plot of log returns and simulated data
	qqplot(stSim[,i], as.numeric(log_returns[,i]), main = paste("Q-Q plot of simulated data with\n", assets[i], "log returns"), pch=20, xlab = paste("Log returns", assets[i]), ylab = "Simulated data");
	abline(a=0,b=1, lty=3)
	q05 = qstable(0.05, stAlpha[i], stBeta[i], stGamma[i], stDelta[i])
	q95 = qstable(0.95, stAlpha[i], stBeta[i], stGamma[i], stDelta[i])
	abline(v=c(q05, 1), lty=2, col='red')
	abline(v=c(q95, 1), lty=2, col='red')
	grid()
}
dev.off()

# # Making two-sample Kolmogorov-Smirnov test
stKSTest = lapply(size, function(p) ks.test(as.numeric(log_returns[,p]), 'pstable', stAlpha[p], stBeta[p], stGamma[p], stDelta[p]))# stSim[,p]))
for(i in size) { 
	stKSTest[[i]]$data.name <- paste("Log returns", assets[i], "and stable distribution");
	stKSTest[[i]]$method = "One-sample Kolmogorov-Smirnov test of stable distribution"; }
stKSTest

# Making Anderson-Darling test
stADtest <- lapply(size, function(p) ad.test(log_returns[,p], pstable, stAlpha[p], stBeta[p], stGamma[p], stDelta[p]))
for(i in size) {
	stADtest[[i]]$data.name <- paste("Log returns", assets[i], "and stable distr");
	stADtest[[i]]$method = "Anderson-Darling GoF test of stable distribution"; }
stADtest

# Perform Crämer-von~Mises test
stCvMtest <- lapply(1:l, function(k) cvm.test(log_returns[,k], 'pstable', stAlpha[k], stBeta[k], stGamma[k], stDelta[k]))
for(i in 1:l) {
	stCvMtest[[i]]$data.name <- paste("Log returns", assets[i])
	stCvMtest[[i]]$method = c("Crämer-von~Mises of GoF", "Null hypothesis: stable distribution") }
stCvMtest

################### MEIXNER DISTRIBUTION ##################

# # Simulate theoretical data
mSim <- sapply(size, function(i) rmeixner(n, mAlpha[i], mBeta[i], mDelta[i], mMu[i])); colnames(mSim)=paste("mSim for",assets)
mSimPDF <- sapply(size, function(i) dmeixner(mSim[,i], mAlpha[i], mBeta[i], mDelta[i], mMu[i]))
colnames(mSimPDF) = paste("mSimPDF for",assets)

# # Show emp histogram with theoretical data line
# pdf(file='Export/Fitting Meixner distribution.pdf', encoding='KOI8-R.enc', width=12, height=7)
setEPS(); postscript("Export/Fitting Meixner distribution.eps", width=12, height=7, paper='special')
for (i in 1:l) {
	par(mfrow = c(1,2))
	# Show emp histogram with theoretical data line
	hist(log_returns[,i], probability = TRUE, label=TRUE, xlab = paste("Log returns", assets[i]), ylim = c(0,1.1*max(mSimPDF[,i])), main = paste(assets[i], "log returns histogram\nwith theoretical PDF"))
	rug(log_returns[,i])
	lines(x, dmeixner(x, mAlpha[i], mBeta[i], mDelta[i], mMu[i]), col="red")
	grid()
	# Show Q-Q plot of log returns and simulated data
	qqplot(mSim[,i], as.numeric(log_returns[,i]), main = paste("Q-Q plot of simulated data with\n", assets[i], "log returns"), pch=20, xlab = paste("Log returns", assets[i]), ylab = "Simulated data");
	abline(a=0,b=1, lty=3)
	q05 = qmeixner(0.05, mAlpha[i], mBeta[i], mDelta[i], mMu[i])
	q95 = qmeixner(0.95, mAlpha[i], mBeta[i], mDelta[i], mMu[i])
	abline(v=c(q05, 1), lty=2, col='red')
	abline(v=c(q95, 1), lty=2, col='red')
	grid()
}
dev.off()

# # Making two-sample Kolmogorov-Smirnov test
mKSTest = lapply(size, function(i) ks.test(as.numeric(log_returns[,i]), 'pmeixner', mAlpha[i], mBeta[i], mDelta[i], mMu[i]))
for(i in size) {
	mKSTest[[i]]$data.name <- paste("Log returns", assets[i], "and Meixner distribution");
	mKSTest[[i]]$method = "One-sample Kolmogorov-Smirnov test of Meixner distribution"; }
mKSTest

# Making Anderson-Darling test
mADtest <- lapply(size, function(p) ad.test(log_returns[,p], pmeixner, mAlpha[p], mBeta[p], mDelta[p], mMu[p]))
for(i in size) {
	mADtest[[i]]$data.name <- paste("Log returns", assets[i], "and Meixner distr");
	mADtest[[i]]$method = "Anderson-Darling GoF test of Meixner distribution"; }
mADtest

# Perform Crämer-von~Mises test
mCvMtest <- lapply(1:l, function(k) cvm.test(log_returns[,k], 'pmeixner', mAlpha[k], mBeta[k], mDelta[k], mMu[k]))
for(i in 1:l) {
	mCvMtest[[i]]$data.name <- paste("Log returns", assets[i])
	mCvMtest[[i]]$method = c("Cramer-von~Mises of GoF", "Null hypothesis: Meixner distribution") }
mCvMtest

######################### SUMMARY #########################

TestResults <- sapply(1:l, function(p) 
  c(hKSTest[[p]]$p.value, stKSTest[[p]]$p.value, mKSTest[[p]]$p.value,
    hADtest[[p]]$p.value, stADtest[[p]]$p.value, mADtest[[p]]$p.value,
    hCvMtest[[p]]$p.value, stCvMtest[[p]]$p.value, mCvMtest[[p]]$p.value))
Dists <- rep(c('H', 'St', 'M'), 3)
Tests <- rep(c('K-S', 'A-D', 'CvM'), each = 3)
rownames(TestResults) <- paste0(Dists, '.', Tests)
colnames(TestResults) <- assets

write.table(round(TestResults,5), 'Export/TestResults.txt', sep='\t', quote=F)
read.table('Export/TestResults.txt', sep='\t')

Marginals <- c('hyperb', 'stable', 'meixner')

best <- sapply(1:l, function(p) which.max(c(hCvMtest[[p]]$p.value, stCvMtest[[p]]$p.value, mCvMtest[[p]]$p.value)))

margins <- Marginals[best]
params <- lapply(1:l, function(i) {switch(margins[i], hyperb = list(Theta = hTheta[,i]), stable = list(alpha=stAlpha[i], beta=stBeta[i], gamma=stGamma[i], delta=stDelta[i]), meixner = list(alpha=mAlpha[i], beta=mBeta[i], delta=mDelta[i], mu=mMu[i]))})


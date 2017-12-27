pdf(file='Export/BestMargins.pdf', width=6, height=12)
# setEPS(); postscript("Export/BestMargins.eps", width=6, height=12)
par(mfrow = c(4,2))
for (i in 1:l) { switch(margins[i],
	hyperb = {
		# Show emp histogram with theoretical data line
		hist(log_returns[,i], probability = TRUE, main = paste(assets[i], "log-returns and\nPDF of the hyperbolic distribution"), xlab = "log-returns", ylim = c(0,1.1*max(hSimPDF[,i])))
		rug(log_returns[,i])
		lines(x, dhyperb(x, hTheta[,i]), col="red")
		grid()
		# Show Q-Q plot of log returns and simulated data
		qqplot(hSim[,i], as.numeric(log_returns[,i]), main=paste("QQ-plot of", assets[i], "log-returns\nand simulated data"), pch=20, xlab = "log-returns", ylab = "Simulated data");
		abline(a=0,b=1, lty=3)
		q05 = qhyperb(0.05, hTheta[,i])
		q95 = qhyperb(0.95, hTheta[,i])
		abline(v=c(q05, 1), lty=2, col='red')
		abline(v=c(q95, 1), lty=2, col='red')
		grid() },
	stable = {
    	# Show emp histogram with theoretical data line
    	hist(log_returns[,i], probability = TRUE, main = paste(assets[i], "log-returns and\nPDF of the stable distribution"), xlab = "log-returns", ylim = c(0,1.1*max(stSimPDF[,i])))
    	rug(log_returns[,i])
    	lines(x, dstable(x, stAlpha[i], stBeta[i], stGamma[i], stDelta[i]), col="red")
    	grid()
    	# Show Q-Q plot of log returns and simulated data
    	qqplot(stSim[,i], as.numeric(log_returns[,i]), main=paste("QQ-plot of", assets[i], "log-returns\nand simulated data"), pch=20, xlab = "log-returns", ylab = "Simulated data");
    	abline(a=0,b=1, lty=3)
    	q05 = qstable(0.05, stAlpha[i], stBeta[i], stGamma[i], stDelta[i])
    	q95 = qstable(0.95, stAlpha[i], stBeta[i], stGamma[i], stDelta[i])
    	abline(v=c(q05, 1), lty=2, col='red')
    	abline(v=c(q95, 1), lty=2, col='red')
    	grid() },
  	meixner = {
    	# Show emp histogram with theoretical data line
    	hist(log_returns[,i], probability = TRUE, main = paste(assets[i], "log-returns and\nPDF of the Meixner distribution"), xlab = "log-returns", ylim = c(0,1.1*max(mSimPDF[,i])))
    	rug(log_returns[,i])
    	lines(x, dmeixner(x, mAlpha[i], mBeta[i], mDelta[i], mMu[i]), col="red")
    	grid()
    	# Show Q-Q plot of log returns and simulated data
    	qqplot(mSim[,i], as.numeric(log_returns[,i]), main=paste("QQ-plot of", assets[i], "log-returns\nand simulated data"), pch=20, xlab = "log-returns", ylab = "Simulated data");
    	abline(a=0,b=1, lty=3)
    	q05 = qmeixner(0.05, mAlpha[i], mBeta[i], mDelta[i], mMu[i])
    	q95 = qmeixner(0.95, mAlpha[i], mBeta[i], mDelta[i], mMu[i])
    	abline(v=c(q05, 1), lty=2, col='red')
    	abline(v=c(q95, 1), lty=2, col='red')
    	grid() }
) }
dev.off()

pdf(file='Export/VaR-ES.pdf', width=11, height=6)
par(mfrow=c(1,2))
plot(level, emp.VaR, xlab=expression(paste(alpha, ' (level)')), ylab='Value-at-Risk', 
     ylim = range(emp.VaR, n.VaR, t.VaR, vc.VaR), main='Value-at-Risk', 'b', lwd=2, lty=3, pch=8); grid()
lines(level, n.VaR, col='blue', 'b', lwd=2, lty=3, pch=24)
lines(level, t.VaR, col='red', 'b', lwd=2, lty=3, pch=25)
lines(level, vc.VaR, col='violet', 'b', lwd=2, lty=3, pch=15)
legend('topleft',col=c('black','blue','red','violet'),legend=c('historical scenario','Gaussian copula', 't-copula','R-vine'), bg='white', lwd=2, lty=3, pch=c(8,24,25))
plot(level, emp.ES, xlab=expression(paste(alpha, ' (level)')), ylab='Expected shortfall', 
     ylim = range(emp.ES, n.ES, t.ES, vc.ES), main='Expected shortfall', 'b', lwd=2, lty=3, pch=8); grid()
lines(level, n.ES, col='blue', 'b', lwd=2, lty=3, pch=24)
lines(level, t.ES, col='red', 'b', lwd=2, lty=3, pch=25)
lines(level, vc.ES, col='violet', 'b', lwd=2, lty=3, pch=15)
legend('topleft',col=c('black','blue','red','violet'),legend=c('historical scenario','Gaussian copula', 't-copula','R-vine'), bg='white', lwd=2, lty=3, pch=c(8,24,25))
dev.off()

pdf(file='Export/VaR-ES-Bootstrap.pdf', width=11, height=6)
# setEPS(); postscript("Export/VaR-ES.eps")
par(mfrow=c(1,2))
plot(level, emp.VaR, main='Value-at-Risk', xlab=expression(paste(alpha, ' (level)')), ylab='Value-at-Risk', ylim = range(emp.VaR, nBootVaR, tBootVaR, boot.VaR), 'b', lty=3, pch=8)
grid()
lines(level, nMeanBootVaR, col='blue', 'b', lty=2, pch=24)
lines(level, tMeanBootVaR, col='red', 'b', lty=2, pch=25)
lines(level, meanBootVaR, col='violet', 'b', lty=2, pch=15)
for (i in 1:length(level)) {
	lines(rep(level[i], 2), lty=3, quantile(nBootVaR[i,], c(0.025, 0.975)), col='blue')
	lines(rep(level[i], 2), lty=3, quantile(tBootVaR[i,], c(0.025, 0.975)), col='red')
	lines(rep(level[i], 2), lty=3, quantile(boot.VaR[i,], c(0.025, 0.975)), col='violet')
}
legend('topleft',col=c('black','blue','red', 'violet'),legend=c('historical scenario','Gaussian copula', "Student's t copula", 'Vine copula'), bg='white', lty=2, pch=c(8,24,25,15))
plot(level, emp.ES, main='Expected Shortfall', xlab=expression(paste(alpha, ' (level)')), ylab='Conditional-Value-at-Risk', ylim = range(emp.ES, nBootES, tBootES, boot.ES), 'b', lty=3, pch=8)
grid()
lines(level, nMeanBootES, col='blue', 'b', lty=2, pch=24)
lines(level, tMeanBootES, col='red', 'b', lty=2, pch=25)
lines(level, meanBootES, col='violet', 'b', lty=2, pch=15)
for (i in 1:length(level)) {
	lines(rep(level[i], 2), lty=3, quantile(nBootES[i,], c(0.025, 0.975)), col='blue')
	lines(rep(level[i], 2), lty=3, quantile(tBootES[i,], c(0.025, 0.975)), col='red')
	lines(rep(level[i], 2), lty=3, quantile(boot.ES[i,], c(0.025, 0.975)), col='violet')
}
legend('topleft',col=c('black','blue','red','violet'),legend=c('historical scenario','Gaussian copula', "Student's t copula", 'Vine copula'), bg='white', lty=2, pch=c(8,24,25,15))
dev.off()

pdf(file='Export/ObservedData.pdf', width=6, height=6)
# setEPS(); postscript("Export/ObservedData.eps")
pairs.panels(log_returns, smooth=F, density=F, ellipses=F, method='kendall', cex.cor=.75)
dev.off()

pdf(file='Export/PseudoObs.pdf', width=6, height=6)
# setEPS(); postscript("Export/PseudoObs.eps")
pairs.panels(empCop, smooth=F, density=F, digits=4, ellipses=F, method='kendall', cex.cor=.75)
dev.off()

# Save results to file
VaRResults <- data.frame(paste0(level*100, '\\%'),
                         paste0('$', 100*round(emp.VaR, 4), '$'),
                         paste0('$', 100*round(n.VaR, 4), '$ & $', 100*round(t.VaR, 4), '$ & $', 100*round(vc.VaR,4), '$'),
                         paste0('$', 1e3*round(nBias.VaR, 5), '$ & $', 1e3*round(tBias.VaR, 5), '$ & $', 1e3*round(vcBias.VaR,5), '$'))

ESResults <- data.frame(paste0(level*100, '\\%'),
                        paste0('$', 100*round(emp.ES, 4), '$'),
                         paste0('$', 100*round(n.ES, 4), '$ & $', 100*round(t.ES, 4), '$ & $', 100*round(vc.ES,4), '$'),
                         paste0('$', 1e3*round(nBias.ES, 5), '$ & $', 1e3*round(tBias.ES, 5), '$ & $', 1e3*round(vcBias.ES,5), '$'))

colnames(VaRResults) <- colnames(ESResults) <- c('Level', 'Empirical, $\\cdot 10^{-2}$', 'Simulated, $\\cdot 10^{-2}$', 'Bias, $\\cdot 10^{-3}$')
write.table(VaRResults, file='Export/VaR Copula Results.txt', sep=' & ', row.names=F, quote=F, eol = " \\\\ \n")
write.table(ESResults, file='Export/ES Copula Results.txt', sep=' & ', row.names=F, quote=F, eol = " \\\\ \n")

# Save results
df1 <- data.frame(paste0(level*100),
                  paste0('$', 100*round(nMeanBootVaR,4), '$ & $', 100*round(tMeanBootVaR,4), '$ & $', 100*round(meanBootVaR,4), '$'),
                  paste0('$', 1e3*round(nBiasBootVaR,5), '$ & $', 1e3*round(tBiasBootVaR,5), '$ & $', 1e3*round(biasBootVaR,5), '$'),
                  paste0('$', 1e3*round(nSDBootVaR,5), '$ & $', 1e3*round(tSDBootVaR,5), '$ & $', 1e3*round(sdBootVaR,5), '$'),
                  paste0('$', 1e3*round(nRMSEBootVaR,5), '$ & $', 1e3*round(tRMSEBootVaR,5), '$ & $', 1e3*round(rmseBootVaR,5), '$'))

df2 <- data.frame(paste0(level*100),
                  paste0('$', 100*round(nMeanBootES,4), '$ & $', 100*round(tMeanBootES,4), '$ & $', 100*round(meanBootES,4), '$'),
                  paste0('$', 1e3*round(nBiasBootES,5), '$ & $', 1e3*round(tBiasBootES,5), '$ & $', 1e3*round(biasBootES,5), '$'),
                  paste0('$', 1e3*round(nSDBootES,5), '$ & $', 1e3*round(tSDBootES,5), '$ & $', 1e3*round(sdBootES,5), '$'),
                  paste0('$', 1e3*round(nRMSEBootES,5), '$ & $', 1e3*round(tRMSEBootES,5), '$ & $', 1e3*round(rmseBootES,5), '$'))

colnames(df1) <- colnames(df2) <- c('Level', 'Mean, $\\cdot 10^{-2}$', 'Bias, $\\cdot 10^{-4}$', 'SD, $\\cdot 10^{-3}$', 'RMSE, $\\cdot 10^{-3}$')

write.table(df1, "Export/VaR. Boot procedure.txt", quote=F, row.names=F, sep=' & ', eol = " \\\\ \n")
write.table(df2, "Export/ES. Boot procedure.txt", quote=F, row.names=F, sep=' & ', eol = " \\\\ \n")
# read.table("Export/VaR. Boot procedure.txt", header=T, sep='&', check.names=F)

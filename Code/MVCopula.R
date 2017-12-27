##############################################
## ---------------- Copula ---------------- ##
##############################################

# ---- Loading files ---- 

source("Testing.r")
source("Functions.r")

# ---- Creating copulas ----

# Fit copula parameters
normFit <- fitCopula(normalCopula(dim=l, dispstr='un'), data=as.matrix(empCop), method='itau')
tFit <- fitCopula(tCopula(dim=l, dispstr='un', df.fixed=FALSE), data=as.matrix(empCop), method='itau')

# Generate the copula
# nCop <- normalCopula(param=Rg, dim=l, dispstr='un')
nCop <- normFit@copula
# tCop <- tCopula(param=Rt, df=dft, dim=l, dispstr='un', df.fixed=T)
tCop <- tFit@copula

# Simulate some points
nSimCop = rCopula(1e3, nCop)
tSimCop = rCopula(1e3, tCop)

# ---- Value-at-Risk and Expected Shortfall ----

# Calculate VaR using copulas and Monte-Carlo
nDataCop <- apply(log_returns, 2, function(p) quantile(p, nSimCop))
n.cop_returns <- nDataCop %*% weights
n.VaR <- quantile(-n.cop_returns, probs=level)

tDataCop <- apply(log_returns, 2, function(p) quantile(p, tSimCop))
t.cop_returns <- tDataCop %*% weights
t.VaR <- quantile(-t.cop_returns, probs=level)

# Bias
nBias.VaR <- n.VaR - emp.VaR
tBias.VaR <- t.VaR - emp.VaR

# Unconditional Kupiec test for VaR
emp.VaRTest <- lapply(1:length(level), function(p) kupiec.test(portRet, emp.VaR[p], level[p]))
n.VaRTest <- lapply(1:length(level), function(p) kupiec.test(n.cop_returns, n.VaR[p], level[p]))
t.VaRTest <- lapply(1:length(level), function(p) kupiec.test(t.cop_returns, t.VaR[p], level[p]))

# Computing Expected Shortfall
n.ES <- sapply(n.VaR, function(p) mean(-n.cop_returns[-n.cop_returns > p]))
t.ES <- sapply(t.VaR, function(p) mean(-t.cop_returns[-t.cop_returns > p]))
# Bias
nBias.ES <- n.ES - emp.ES
tBias.ES <- t.ES - emp.ES

# # # # # # # # # # # # # # # # # # # # # # # # 
# Create new file to save plots
pdf(file="Export/Multivariate Gaussian and Student's t copula.pdf", encoding='KOI8-R.enc', width=12, height=7)
# setEPS(); postscript("Export/Multivariate Gaussian and Student's t copula.eps")
# Plot VaR and ES
par(mfrow=c(1,2))
plot(level, emp.VaR, xlab=expression(paste(alpha, ' (level)')), ylab='Value-at-Risk', 
     ylim = range(emp.VaR, n.VaR, t.VaR), main='Value-at-Risk', 'b', lwd=2, lty=3, pch=8); grid()
lines(level, n.VaR, col='blue', 'b', lwd=2, lty=3, pch=24)
lines(level, t.VaR, col='red', 'b', lwd=2, lty=3, pch=25)
legend('topleft',col=c('black','blue','red'),legend=c('historical scenario','Gaussian copula', 't-copula'), bg='white', lty=3, pch=c(8,24,25))
plot(level, emp.ES, xlab=expression(paste(alpha, ' (level)')), ylab='Expected shortfall', 
     ylim = range(emp.ES, n.ES, t.ES), main='Expected shortfall', 'b', lwd=2, lty=3, pch=8); grid()
lines(level, n.ES, col='blue', 'b', lwd=2, lty=3, pch=24)
lines(level, t.ES, col='red', 'b', lwd=2, lty=3, pch=25)
legend('topleft',col=c('black','blue','red'),legend=c('historical scenario','Gaussian copula', 't-copula'), bg='white', lty=3, pch=c(8,24,25))
dev.off()
# # # # # # # # # # # # # # # # # # # # # # # # 

# ---- Bootstrap VaR and ES ----

# Genearate copula samples
# nCopulaBoot <- boot.nCopula(log_returns)
nCopulaBoot <- replicate(200, rCopula(n, nCop))
# tCopulaBoot <- boot.tCopula(log_returns)
tCopulaBoot <- replicate(200, rCopula(n, tCop))

# Calculate VaR using bootstrap copulas and Monte-Carlo
nBootVaR <- VaR.bootstrap(log_returns, nCopulaBoot, weights, level)
tBootVaR <- VaR.bootstrap(log_returns, tCopulaBoot, weights, level)

# Calculate mean, bias, sd, mse
nMeanBootVaR <- apply(nBootVaR, 1, function(p) mean(p))
tMeanBootVaR <- apply(tBootVaR, 1, function(p) mean(p))
nDiffBootVaR <- sapply(1:length(probs), function(p) max(abs(nMeanBootVaR[p] - nBootVaR[p,])))
tDiffBootVaR <- sapply(1:length(probs), function(p) max(abs(tMeanBootVaR[p] - tBootVaR[p,])))

nBiasBootVaR <- sapply(1:length(probs), function(p) nMeanBootVaR[p] - emp.VaR[p])
tBiasBootVaR <- sapply(1:length(probs), function(p) tMeanBootVaR[p] - emp.VaR[p])
nSDBootVaR <- apply(nBootVaR, 1, function(p) sd(p))
tSDBootVaR <- apply(tBootVaR, 1, function(p) sd(p))
nRMSEBootVaR <- sapply(1:length(probs), function(p) sqrt(1/ncol(nBootVaR)*sum((nBootVaR[p,] - emp.VaR[p])^2)))
tRMSEBootVaR <- sapply(1:length(probs), function(p) sqrt(1/ncol(tBootVaR)*sum((tBootVaR[p,] - emp.VaR[p])^2)))

# Calculate ES for each pobs array
nBootES <- ES.bootstrap(log_returns, nCopulaBoot, weights, level)
tBootES <- ES.bootstrap(log_returns, tCopulaBoot, weights, level)

# Calculate mean, bias, sd, mse
nMeanBootES <- apply(nBootES, 1, function(p) mean(p))
tMeanBootES <- apply(tBootES, 1, function(p) mean(p))
nDiffBootES <- sapply(1:length(probs), function(p) max(abs(nMeanBootES[p] - nBootES[p,])))
tDiffBootES <- sapply(1:length(probs), function(p) max(abs(tMeanBootES[p] - tBootES[p,])))

nBiasBootES <- sapply(1:length(probs), function(p) nMeanBootES[p]-emp.ES[p])
tBiasBootES <- sapply(1:length(probs), function(p) tMeanBootES[p]-emp.ES[p])
nSDBootES <- apply(nBootES, 1, function(p) sd(p))
tSDBootES <- apply(tBootES, 1, function(p) sd(p))
nRMSEBootES <- sapply(1:length(probs), function(p) sqrt(1/ncol(nBootES)*sum((nBootES[p,] - emp.ES[p])^2)))
tRMSEBootES <- sapply(1:length(probs), function(p) sqrt(1/ncol(tBootES)*sum((tBootES[p,] - emp.ES[p])^2)))

# # # # # # # # # # # # # # # # # # # # # # # # 
# Plot VaR and ES
pdf(file='Export/Gaussian and t-copula. Boot VaR & ES.pdf', encoding='KOI8-R.enc', width=12, height=7)
par(mfrow=c(1,2))
plot(level, emp.VaR, xlab=expression(paste(alpha, ' (level)')), ylab='Value-at-Risk', 
     ylim = range(emp.VaR, nBootVaR, tBootVaR), main='Value-at-Risk\nBootstrap procedure', 
     'b', lwd=2, lty=3, pch=8)
grid()
lines(level, nMeanBootVaR, col='blue', 'b', lwd=2, lty=3, pch=24)
lines(level, tMeanBootVaR, col='red', 'b', lwd=2, lty=3, pch=25)
for (i in 1:length(level)) {
	lines(rep(level[i], 2), quantile(nBootVaR[i,], c(0.025, 0.975)), lwd=2, col='blue')
	lines(rep(level[i], 2), quantile(tBootVaR[i,], c(0.025, 0.975)), lwd=2, col='red')
}
legend('topleft', legend=c('historical scenario','Gaussian copula', 't-copula'), 
       col=c('black','blue','red'), bg='white', lwd=2, lty=3, pch=c(8,24,25))
plot(level, emp.ES, xlab=expression(paste(alpha, ' (level)')), ylab='Expected shortfall', 
     ylim = range(emp.ES, nBootES, tBootES), main='Expected shortfall\nBootstrap procedure', 
     'b', lwd=2, lty=3, pch=8)
grid()
lines(level, nMeanBootES, col='blue', 'b', lwd=2, lty=3, pch=24)
lines(level, tMeanBootES, col='red', 'b', lwd=2, lty=3, pch=25)
for (i in 1:length(level)) {
	lines(rep(level[i], 2), quantile(nBootES[i,], c(0.025, 0.975)), lwd=2, col='blue')
	lines(rep(level[i], 2), quantile(tBootES[i,], c(0.025, 0.975)), lwd=2, col='red')
}
legend('topleft', legend=c('historical scenario','Gaussian copula', 't-copula'),
       col=c('black','blue','red'), bg='white', lwd=2, lty=3, pch=c(8,24,25))
dev.off()
# # # # # # # # # # # # # # # # # # # # # # # # 

# Copula goodness-of-fit test
nCopGoF <- gofCopula(nCop, as.matrix(log_returns), N=100)
tCopGoF <- gofCopula(tCop, as.matrix(log_returns), N=100)


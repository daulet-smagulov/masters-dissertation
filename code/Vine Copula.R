#################################################
## --------------- Vine Copula --------------- ##
#################################################

# ---- Loading libraries ---- 

if(!require(devtools)) {install.packages("devtools"); require(devtools)}
if(!require(VineCopula)) {devtools::install_github("tnagler/VineCopula"); require(VineCopula)}
if(!require(pracma)) {install.packages("pracma"); require(pracma)}
if(!require(psych)) {install.packages("psych"); require(psych)}

source("Testing.r")
source("Functions.r")

# ---- Create copula ----

# Type of vine
type = "RVine"

# Structure of copula
Tree <- RVineStructureSelect(data=empCop, type=type, treecrit = 'rho')

# Generate vine copula and multivariate distribution objects
Copula <- vineCopula(Tree)
Mvdc <- mvdc(copula=Copula, margins, paramMargins=params)

# GofTest <- RVineGofTest(empCop, Tree)

# Simulate copula and scenarios
R <- scenarioGenerate(1000, Mvdc)
simCop <- RVineSim(1e3, Tree)
# simCop <- rCopula(1e3, Copula)
# simCop <- pobs(R)

# # # # # # # # # # # # # # # # # # # # # # # # 
# Save file
pdf(file='Export/Vine Copula.pdf', encoding='KOI8-R.enc', width=8.5, height=9)
# setEPS(); postscript("Export/Vine Copula.eps")
for (i in 1:(l-1)) plot(Tree, tree=i)
pairs.panels(empCop, smooth=F, density=F, digits=4, ellipses=F, method='kendall', main='Empirical Copula', cex.cor=.75)
pairs.panels(simCop, smooth=F, density=F, digits=4, ellipses=F, method='kendall', main='Simulated Copula', cex.cor=.75)
pairs.panels(log_returns, smooth=F, density=F, ellipses=F, method='kendall', main='Observed data', cex.cor=.75)
pairs.panels(R, smooth=F, density=F, ellipses=F, method='kendall', main='Simulated data', cex.cor=.75)
dev.off()
# # # # # # # # # # # # # # # # # # # # # # # # 

# ---- VaR and CVaR ----

# Calculate VaR using copulas and Monte-Carlo
DataCop <- apply(log_returns, 2, function(p) quantile(p, simCop))
cop_returns <- DataCop %*% weights
vc.VaR <- quantile(-cop_returns, probs=level)

# Bias
vcBias.VaR <- vc.VaR - emp.VaR

# Kupiec test of VaR
empVaRTest <- lapply(1:length(level), function(p) kupiec.test(portRet, emp.VaR[p], level[p]))
vcVaRTest <- lapply(1:length(level), function(p) kupiec.test(cop_returns, vc.VaR[p], level[p]))

# Computing Expected Shortfall
vc.ES <- sapply(vc.VaR, function(p) mean(-cop_returns[-cop_returns > p]))

# Bias
vcBias.ES <- vc.ES - emp.ES

# # # # # # # # # # # # # # # # # # # # # # # # 
# Create the file
pdf(file='Export/Vine Copula. VaR & ES.pdf', encoding='KOI8-R.enc', width=12, height=7)
# setEPS(); postscript("Export/Vine Copula. VaR & ES.eps")
# Plot VaR and ES
par(mfrow=c(1,2))
plot(level, emp.VaR, xlab=expression(paste(alpha, ' (level)')), ylab='Value-at-Risk', ylim = range(emp.VaR, vc.VaR), main='Value-at-Risk', 'b', lwd=2, lty=3, pch=8)
grid()
lines(level, vc.VaR, col='violet', 'b', lwd=2, lty=3, pch=24)
legend('bottomright', col=c('black', 'violet'), legend=c('historical scenario','R-vine copula'), bg='white', lwd=2, lty=3, pch=c(8,24))
plot(level, emp.ES, xlab=expression(paste(alpha, ' (level)')), ylab='Expected shortfall', ylim = range(emp.ES, vc.ES), main='Expected shortfall', 'b', lwd=2, lty=3, pch=8)
grid()
lines(level, vc.ES, col='violet', 'b', lwd=2, lty=3, pch=24)
legend('bottomright', col=c('black', 'violet'), legend=c('historical scenario','R-vine copula'), bg='white', lwd=2, lty=3, pch=c(8,24))
# # # # # # # # # # # # # # # # # # # # # # # # 

# ---- Bootstrap VaR and ES ----

# Genearate copula samples
# bootVineCop <- bootVineCopula(log_returns)
bootVineCop <- replicate(200, RVineSim(n, Tree))

# Calculate VaR for each pobs array
boot.VaR <- VaR.bootstrap(log_returns, bootVineCop, weights, level)

# Calculate mean, bias, sd, mse
meanBootVaR <- apply(boot.VaR, 1, function(p) mean(p))

biasBootVaR <- sapply(1:length(probs), function(p) meanBootVaR[p]-emp.VaR[p])
sdBootVaR <- apply(boot.VaR, 1, function(p) sd(p))
rmseBootVaR <- sapply(1:length(probs), function(p) sqrt(1/ncol(boot.VaR)*sum((boot.VaR[p,] - emp.VaR[p])^2)))

# Calculate ES for each pobs array
boot.ES <- ES.bootstrap(log_returns, bootVineCop, weights, level)

# Calculate mean, bias, sd, mse
meanBootES <- apply(boot.ES, 1, function(p) mean(p))

biasBootES <- sapply(1:length(probs), function(p) meanBootES[p]-emp.ES[p])
sdBootES <- apply(boot.ES, 1, function(p) sd(p))
rmseBootES <- sapply(1:length(probs), function(p) sqrt(1/ncol(boot.ES)*sum((boot.ES[p,] - emp.ES[p])^2)))

# # # # # # # # # # # # # # # # # # # # # # # # 
# Plot VaR and ES
# pdf(file='Export/Boot.VaR.ES.pdf', encoding='KOI8-R.enc', width=12, height=7)
par(mfrow=c(1,2))
plot(level, emp.VaR, xlab=expression(paste(alpha, ' (level)')), ylab='Value-at-Risk', ylim = range(emp.VaR, boot.VaR), main='Value-at-Risk\nBootstrap procedure', 'b', lwd=2, lty=3, pch=8)
grid()
lines(level, meanBootVaR, col='violet', 'b', lwd=2, lty=3, pch=15)
for (i in 1:length(level)) 
lines(rep(level[i], 2), quantile(boot.VaR[i,], c(0.025, 0.975)), lwd=2, col='violet')
legend('bottomright',col=c('black', 'violet'),legend=c('historical scenario','R-vine copula'), bg='white', lwd=2, lty=3, pch=c(8,15))
plot(level, emp.ES, xlab=expression(paste(alpha, ' (level)')), ylab='Expected shortfall', ylim = range(emp.ES, boot.ES), main='Expected shortfall\nBootstrap procedure', 'b', lwd=2, lty=3, pch=8)
grid()
lines(level, meanBootES, col='violet', 'b', lwd=2, lty=3, pch=15)
for (i in 1:length(level)) 
lines(rep(level[i], 2), quantile(boot.ES[i,], c(0.025, 0.975)), lwd=2, col='violet')
legend('bottomright',col=c('black', 'violet'),legend=c('historical scenario','R-vine copula'), bg='white', lwd=2, lty=3, pch=c(8,15))
dev.off()
# # # # # # # # # # # # # # # # # # # # # # # # 


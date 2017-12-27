# rm(list=ls(all=TRUE))
# setwd("/Users/Daulet/Documents/University/Masters/Science/SP/Main")
if(!require(parallel)) {install.packages("parallel"); require(parallel)}
# load the data from "Finam"
if(!require(rusquant)) {install.packages("rusquant", repos="http://R-Forge.R-project.org", type="source"); library(rusquant)}
if(!require(curl)) {install.packages("curl"); library(curl)}
# Distributions
if(!require(HyperbolicDist)) {install.packages("HyperbolicDist"); library(HyperbolicDist)}
if(!require(stabledist)) {install.packages("stabledist"); library(stabledist)}
if(!require(Runuran)) {install.packages("Runuran"); library(Runuran)}
# Fit distr by PDF and CDF
if(!require(fitdistrplus)) {install.packages("fitdistrplus"); library(fitdistrplus)}
# For parameters esttimation via "Method of moments"
if(!require(moments)) {install.packages("moments"); library(moments)}
# Portfolio optimization
if(!require(PortfolioAnalytics)) {install.packages("PortfolioAnalytics"); library(PortfolioAnalytics)}
if(!require(ROI)) {install.packages("ROI"); library(ROI)}
if(!require(ROI.plugin.glpk)) {install.packages("ROI.plugin.glpk"); library(ROI.plugin.glpk)}
if(!require(PerformanceAnalytics)) {install.packages("PerformanceAnalytics"); library(PerformanceAnalytics)}
# Copula main package
if(!require(copula)) {install.packages("copula"); library(copula)} 
# if(!require(parallel)) {install.packages("parallel"); library(parallel)} 

source("Functions.r")

###################### PREPARING DATA #####################

set.seed(20150509) # rm(.Random.seed, envir=globalenv())
assets <- c("RTS", "SBRF", "GAZR", "GMKR")

# Data for the last year
# start <- Sys.Date() - 365; end <- Sys.Date()
start <- '2015-12-16'; end <- '2017-12-16'
getSymbols(paste0('SPFB.', assets), src='Finam', from=start, to=end, auto.assign=TRUE)

price = merge.zoo(SPFB.RTS[,4], SPFB.SBRF[,4], SPFB.GAZR[,4], SPFB.GMKR[,4])

##################### INITIALIZATION ######################

log_returns <- diff(log(price))
colnames(log_returns) <- assets
date <- as.Date(rownames(as.matrix(log_returns)))

n <- nrow(log_returns)
l <- length(assets)
size <- seq(1, l)

Dims <- c(max(floor(l/2),1), ceiling(l/max(floor(l/2),1)))
dims <- c(min(Dims), max(Dims))

# # Draw histograms
Hist = lapply(size, function(p) hist(log_returns[,p], plot = FALSE))

pdf(file='Export/Histogram.pdf', encoding='KOI8-R.enc', width=10, height=12)
# setEPS(); postscript("Export/Histogram.eps")
# quartz("Log returns histogram"); par(mfrow = dims);
# windows(title="Hyperbolic distribution"); 
par(mfrow = dims);
for (p in size) {
    # pdf(file=paste0('Export/Histogram_',p,'.pdf'), encoding='KOI8-R.enc', width=8, height=9)
	hist(log_returns[,p], probability = TRUE, label=TRUE, xlab = paste("Log returns", assets[p]), ylim=c(0,1.1*max(Hist[[p]]$density)), main = paste(assets[p])); 
	rug(log_returns[,p]);
	grid()
	# dev.off()
}
dev.off()

# Finding moments
meanLR = apply(log_returns, 2, mean)#; names(meanLR) = assets
sdLR = apply(log_returns, 2, sd)
varLR = sapply(size, function(i) var(log_returns[,i]))#; names(varLR) = assets
skwLR = sapply(size, function(i) skewness(log_returns[,i]))#; names(skwLR) = assets
krtLR = sapply(size, function(i) 3+kurtosis(log_returns[,i]))#; names(krtLR) = assets
medLR = sapply(size, function(i) quantile(log_returns[,i], 0.5))#; names(medLR) = assets

spearmans_rho = cor(log_returns, method='spearman')
kendalls_tau = cor(log_returns, method='kendall')

# Save moments
write.table(round(cbind(meanLR, sdLR, spearmans_rho, kendalls_tau), 3), file='Export/Assets.txt', quote=FALSE)

# Parallel preprocessing ----
# Calculate number of cores
# no_cores <- detectCores() - 1

##################### FIT THE PARAMETERS ######################

# ---- Hyperbolic distribution ----

# # Fit a hyperbolic parameter distribution to data and display the results
hFit <- lapply(size, function(p) hyperbFit(log_returns[,p]))
for(i in size) hFit[[i]]$obsName <- paste("Hyperbolic dist parameters for", assets[i], "log returns")

# # Finding fitted Theta
hTheta <- sapply(size, function(p) hFit[[p]]$Theta)
colnames(hTheta) = assets

# # Calculate standard errors for the estimates of pi, zeta, delta, and mu
# summary(hyperbFit(log_returns[,1],hessian=TRUE))
# summary(hyperbFit(log_returns[,2],hessian=TRUE))

# ---- Stable distribution ----

# library(fBasics) # fit the parameters' estimates (the fastest one)

# # Making dist and prob functions for fitting
# # Note that
# alpha = 2*exp(a)/(1+exp(a))
# beta = (exp(b)-1)/(exp(b)+1)
# gamma = exp(c)
# delta = d

# Fit the parameters for maximum likelihood

# # Using "fitdistrplus" package, maximum gooness-of-fit estimation, 
# # Cramer-von Mises distance
# stTheta = lapply(size, function(i) fitdist(as.numeric(log_returns[,i]), "stfit", list(a=1, b=0, c=1, d=0), method = "mge", gof="CvM"))
stTheta = apply(log_returns, 2, fitdist, "stfit", list(a=1, b=0, c=1, d=0), 
                method = "mge", gof="CvM")

# # Initiate cluster
# cl <- makeCluster(no_cores)
# clusterExport(cl, list('fitdist', 'dstfit', 'pstfit', 'qstfit', 'log_returns'))
# stTheta = parApply(cl, log_returns, 2, fitdist, "stfit",
#                    list(a=1, b=0, c=1, d=0),
#                    method="mge", gof="CvM")
# stopCluster(cl)
# # summary(stTheta[[1]])
# # stTheta[[2]]$estimate

# # Parameters' estimates
stAlpha = sapply(size, function(i) 2*exp(stTheta[[i]]$estimate[1])/(1+exp(stTheta[[i]]$estimate[1]))); names(stAlpha) = assets
stBeta = sapply(size, function(i) 2*exp(stTheta[[i]]$estimate[2])/(1+exp(stTheta[[i]]$estimate[2]))-1); names(stBeta) = assets
stGamma = sapply(size, function(i) exp(stTheta[[i]]$estimate[3])); names(stGamma) = assets
stDelta = sapply(size, function(i) stTheta[[i]]$estimate[4]); names(stDelta) = assets

# Parameter' estimates:
# stAlpha <- c(1.86195475042, 1.80090675045, 1.85289686513, 1.75151042158)
# stBeta <- c(0.193447086091, 0.909624863815, 0.843236477691, 0.873821610727)
# stGamma <- c(0.0123612530950, 0.0101918408078, 0.0081610302912, 0.0092034946389)
# stDelta <- c(0.001595556, 0.0007162329, -0.0005917765, -0.001533033)

# stAlpha <- c(1.86809613117, 1.88720434447, 1.91313305721, 1.99999382891)
# stBeta <- c(.779109526463, .104567242325, .893104350557, -.949878849627)
# stGamma <- c(.00942246796194, .00983968276316, .00723586854282, .00952363946488)
# stDelta <- c(.00023414912761, .00106219313656,-.00144839630822,-.00060076760893)

names(stAlpha) <- names(stBeta) <- names(stGamma) <- names(stDelta) <- assets

# ---- Meixner distribution ----

# # Finding parameters' estimates by empirical moments
mAlpha.start = sqrt(varLR*(2*krtLR-3*skwLR^2-6))
mBeta.start = sign(skwLR)*acos((krtLR-2*skwLR^2-3)/(krtLR-skwLR^2-3))
mDelta.start = 1/(krtLR-skwLR^2-3)
mMu.start = as.numeric(meanLR-mAlpha.start*mDelta.start*tan(mBeta.start/2))

# # # Finding parameters' estimates by "fitdist" function
mTheta = lapply(size, function(k) fitdist(as.numeric(log_returns[,k]), "meixner", 
                                          list(alpha=mAlpha.start[k], beta=mBeta.start[k], 
                                               delta=mDelta.start[k], mu=mMu.start[k]), 
                                          lower = c(0, -pi, 0, -Inf), upper=c(2,pi,Inf,Inf), 
                                          method = "mge", gof="CvM"))

# Initiate cluster
# cl <- makeCluster(no_cores)
# mTheta = parLapply(cl, size, function(k) fitdist(as.numeric(log_returns[,k]), "meixner",
#                                                  list(alpha=mAlpha.start[k], beta=mBeta.start[k], 
#                                                       delta=mDelta.start[k], mu=mMu.start[k]), 
#                                                  lower = c(0, -pi, 0, -Inf), upper=c(2,pi,Inf,Inf), 
#                                                  method = "mge", gof="CvM"))
# stopCluster(cl)

mAlpha = sapply(size, function(p) mTheta[[p]]$estimate[1])
mBeta = sapply(size, function(p) mTheta[[p]]$estimate[2])
mDelta = sapply(size, function(p) mTheta[[p]]$estimate[3])
mMu = sapply(size, function(p) mTheta[[p]]$estimate[4])

# mAlpha <- c(.0076304543206, .0215295303468, .0245375426974, .0016884968496)
# mBeta <- c(.99240464777, .36726529709, .08161891023, 1.43730618121)
# mDelta <- c(5.0371543521, .9750308539, .4813272349, 72.2165969972)
# mMu <- c(-.019338674255, -.002347532776, -.001393296757, -.107099677171)

# mAlpha <- c(0.0192546766828,   0.0272593124369,  0.0220921896274,  0.0046423569107)
# mBeta  <- c(0.4721389574,      0.6822235316,     0.2667767709,     2.1214394827)
# mDelta <- c(1.7461305524,      0.6724510759,     0.6980944439,     4.2282651257)
# mMu    <- c(-0.0057512291368, -0.0034474201586, -0.0015772668879, -0.0341261894181)

names(mAlpha) <- names(mBeta) <- names(mDelta) <- names(mMu) <- assets

# mAlpha = sapply(size, function(p) sqrt(varLR[p]*(2*krtLR[p]-3*skwLR[p]^2-6)))
# mBeta = sapply(size, function(p) sign(skwLR[p])*acos((krtLR[p]-2*skwLR[p]^2-3)/(krtLR[p]-skwLR[p]^2-3)))
# mDelta = sapply(size, function(p) 1/(krtLR[p]-skwLR[p]^2-3))
# mMu = sapply(size, function(p) meanLR[p]-mAlpha1[p]*mDelta1[p]*tan(mBeta1[p]/2))

# ---- Save parameters ----

StTheta <- rbind(stAlpha, stBeta, stGamma, stDelta)
MTheta <- rbind(mAlpha, mBeta, mDelta, mMu)

Parameters <- rbind(hTheta, StTheta, MTheta)


write.table(round(Parameters, 5), 'Export/Parameters.txt', sep='\t', quote=F)
read.table('Export/Parameters.txt', sep='\t')

# # # # # # # # # # # # # # # # # # # # # # # # 

# Looking for optimal portfolio
# level <- c(.999, .995, .99, .95, .9) #; level <- 1-probs
level <- c(.9, .95, .99, .995, .999)

# Create new portfolio with l assets
portfolio <- portfolio.spec(assets)
# Add "box" and "full-investment" constraints
portfolio <- add.constraint(portfolio, type='box', min=.05, max=1)
# portfolio <- add.constraint(portfolio, type='long_only')
portfolio <- add.constraint(portfolio, type='full_investment')
# Minimize ES
min.ES <- add.objective(portfolio, type='risk', name='ES')
# Optimize the portfolio
opt.ES <- optimize.portfolio(log_returns, portfolio=min.ES, optimize_method='ROI', trace=T)
weights <- opt.ES$weights

# ---- Empirical Value-at-Risk and Expected Shortfall ----

# Calculate empirical VaR and ES
portRet <- Return.portfolio(log_returns, weights)
emp.VaR <- quantile(-portRet, probs=level)
emp.ES <- sapply(emp.VaR, function(p) mean(-portRet[-portRet > p]))

eq.weights <- rep(1/l, l)
eq.portRet <- Return.portfolio(log_returns, eq.weights)
eq.emp.VaR <- quantile(-eq.portRet, probs=level)
eq.emp.ES <- sapply(eq.emp.VaR, function(p) mean(-eq.portRet[-eq.portRet > p]))

eq.bias.VaR <- eq.emp.VaR - emp.VaR 
eq.bias.ES <- eq.emp.ES - emp.ES

VaR.ES <- data.frame(paste0(level*100, '\\%'), 
	paste0('$', 100*round(emp.VaR, 4), '$ & $', 100*round(emp.ES, 4), '$'),
    paste0('$', 100*round(eq.emp.VaR, 4), '$ & $', 100*round(eq.emp.ES, 4), '$'),
    paste0('$', 100*round(eq.bias.VaR, 4), '$ & $', 100*round(eq.bias.ES, 4), '$'))
colnames(VaR.ES) <- c('Level', 'Optimal portfolio VaR / ES, $\\cdot 10^{-2}$', 'Equiweighted portfolio VaR / ES, $\\cdot 10^{-2}$', 'Bias, $\\cdot 10^{-2}$')
write.table(VaR.ES, "Export/VaR.ES.txt", row.names = F, quote = F, sep = " & ", eol = " \\\\ \n")


# ---- Empirical copula ----

empCop <- pobs(log_returns)


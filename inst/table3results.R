# This code produces the results in Table 3 of the manuscript:
# "Transmission thresholds for the spread of infections in healthcare facilities"

rm(list=ls())

snstvty <- 0.85
gam <- 1/387
eps <- 0.5

mu <- 33.8
l25 <- 16
l50 <- 28
l75 <- 43

z95 <- qnorm(0.975)

getNorm <- function(n,ns) {p <- ns/n; c(mean=p, sd=sqrt(p*(1-p)/n))}

haydenSet <- rbind(admPos = c(mean=0.206, sd=0.017/z95),
                   xSecPos = c(mean=0.458, sd=0.037/z95),
                   clinDetInc = getNorm(n=178516, ns=656))

rLHS <- read.table("inst/extdata/rLHS.txt", header=TRUE)
numRuns <- nrow(rLHS)

getLHSval <- function(prmName,i)
	as.numeric(qnorm(rLHS[i,prmName], mean=haydenSet[prmName,'mean'], sd=haydenSet[prmName,'sd']))

cdf <- function(tm,px,rx,rg,k) ifelse(rx>0 & rg>0, px*pexp(tm,rx) + (1-px)*pgamma(tm,shape=k,rate=rg), 0)
getMu <- function(px,rx,rg,k) px/rx + (1-px)*k/rg

dat  <- c(mu,l25,l50,l75)

losErrFn <- function(mu,l25,l50,l75,px,rx,rg,k){
	muEst <- getMu(px,rx,rg,k)
	p25 <- cdf(l25,px,rx,rg,k)
	p50 <- cdf(l50,px,rx,rg,k)
	p75 <- cdf(l75,px,rx,rg,k)

	sum((c(muEst,100*p25,100*p50,100*p75)-c(mu,25,50,75))^2)
}

parInit <- c(px = 0.5, rx = 1/40, rg = 5/30, k = 5)
losSol <- optim(fn=function(x) do.call(losErrFn,as.list(c(dat,x))), par=parInit, control=list(reltol=1e-18,maxit=10000))

losPrm <- losSol$par

MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
	sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

mgf <- function(x, deriv=0) MGFmixedgamma(x, prob = c(losPrm['px'],1-losPrm['px']),
                                          rate = c(losPrm['rx'],losPrm['rg']),
                                          shape = c(1,losPrm['k']), deriv)

Sfun <- function(a) -a
Afun <- function(a) rbind(a,0)
Cfun <- function(d) rbind(c(-gam-d,0),c(d,0))
R <- cbind(gam,0)
transm <- c(1,1-eps)
colStates <- 2:3
preClinStates <- 2

getCalibModel3 <- function(p){

	init <- c(1-p['PA'],p['PA'],0)

	getequilib <- function(statePrm){
	  C <- Cfun(statePrm['deltaC'])
	  S <- Sfun(statePrm['alpha'])
	  A <- Afun(statePrm['alpha'])
	  W <- rbind(cbind(S,R),cbind(A,C))
	  equilib(W,init,mgf)
	}

	errFn <- function(statePrm){
		eq <- getequilib(statePrm)
		eqPrev <- sum(eq[colStates]) * 100
		eqClin <- sum(eq[preClinStates]) * statePrm['deltaC'] * 10000
		sum((c(eqPrev,eqClin) - c(p['P']*100, p['clinDetInc']*10000))^2)
	}

	sol <- optim(c(alpha = 0.018, deltaC = 0.008), errFn, control=list(reltol=1e-18,maxit=10000))
	solState <- sol$par
	solEq <- getequilib(solState)
	solBeta <- as.numeric(solState['alpha']/sum(solEq[colStates]*transm))
	R0 <- facilityR0(S = Sfun(0),
	                 C = Cfun(solState['deltaC']),
	                 A = Afun(1),
	                 transm = solBeta*transm,
	                 initS = 1,
	                 mgf = mgf)

	c(beta=solBeta,solState['deltaC'],R0=as.numeric(R0))
}

get95ci <- function(d){
	i1 <- round(0.025*numRuns) + 1
	i2 <- numRuns - round(0.025*numRuns)
	apply(d,2,sort)[c(i1,i2),]
}

p <- c(PA = haydenSet['admPos','mean'] / snstvty,
       P = haydenSet['xSecPos','mean'] / snstvty,
       clinDetInc = haydenSet['clinDetInc','mean'])

meanModel3 <- getCalibModel3(p)

meanTable <- matrix(NA,7,1)
rownames(meanTable) <- c(names(losPrm),names(meanModel3))
colnames(meanTable) <- 'Estimate'
meanTable[names(losPrm),1] <- losPrm
meanTable[names(meanModel3),1] <- meanModel3

mcModel3 <- matrix(0,numRuns,3)
colnames(mcModel3) <- names(meanModel3)

for(n in 1:numRuns){
	p <- c(PA = getLHSval('admPos',n) / snstvty,
	       P = getLHSval('xSecPos',n) / snstvty,
	       clinDetInc = getLHSval('clinDetInc',n))
	mcModel3[n,] <- getCalibModel3(p)
}

ciModel3 <- get95ci(mcModel3)


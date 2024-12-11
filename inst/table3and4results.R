# This code produces the results in Tables 3 and 4 of the manuscript:
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

table3mean <- matrix(NA,7,1)
rownames(table3mean) <- c(names(losPrm),names(meanModel3))
colnames(table3mean) <- 'Estimate'
table3mean[names(losPrm),1] <- losPrm
table3mean[names(meanModel3),1] <- meanModel3

print(table3mean)

mcModel3 <- matrix(0,numRuns,3)
colnames(mcModel3) <- names(meanModel3)

for(n in 1:numRuns){
	p <- c(PA = getLHSval('admPos',n) / snstvty,
	       P = getLHSval('xSecPos',n) / snstvty,
	       clinDetInc = getLHSval('clinDetInc',n))
	mcModel3[n,] <- getCalibModel3(p)
}

table3CI <- t(get95ci(mcModel3))
colnames(table3CI) <- c("low95CI","high95CI")
print(table3CI)

getR0intervention <- function(bet, dc, ds, gamd){
  if(gamd < Inf){
    return(facilityR0(S = diag(0,2),
                      C = rbind(c(-ds-dc-gam,0,0),
                                c(ds,-dc-gamd,0),
                                c(dc,dc,0)),
                      A = rbind(c(1,0),c(0,1-eps),c(0,0)),
                      transm = bet*c(1,1-eps,1-eps),
                      initS = rbind(1,0), mgf))
  }
  facilityR0(S = diag(0,2),
             C = rbind(c(-ds-dc-gam,0),c(dc,0)),
             A = rbind(c(1,0),c(0,0)),
             transm = bet*c(1,1-eps),
             initS = rbind(1,0), mgf)
}
getR0intervention <- Vectorize(getR0intervention)

surveillancethreshold <- function(bet, dc, gamd){

  if(getR0intervention(bet=bet,dc=dc,ds=0,gamd=gamd) < 1) return(0)

  fn <- function(ds) (getR0intervention(bet, dc, ds, gamd) - 1)^2
  optimize(fn, c(0,1), tol=1e-10)$minimum
}
surveillancethreshold <- Vectorize(surveillancethreshold)

dsweekly <- 1/7 * 0.954 * 0.85
dsbiweekly <- 1/14 * 0.954 * 0.85

getTable4row <- function(gamd){

  gamdtable <- c(decolratefactor = gamd/gam, meandaystodecol = 1/gamd)
  meanR0weekly <- getR0intervention(bet = as.numeric(meanModel3["beta"]),
                                    dc = as.numeric(meanModel3["deltaC"]),
                                    ds = dsweekly,
                                    gamd = gamd)

  mcR0weekly <- getR0intervention(bet = as.numeric(mcModel3[,"beta"]),
                                  dc = as.numeric(mcModel3[,"deltaC"]),
                                  ds = dsweekly,
                                  gamd = gamd)

  ciR0weekly <- get95ci(cbind(mcR0weekly))

  meanR0biweekly <- getR0intervention(bet = as.numeric(meanModel3["beta"]),
                                      dc = as.numeric(meanModel3["deltaC"]),
                                      ds = dsbiweekly,
                                      gamd = gamd)

  mcR0biweekly <- getR0intervention(bet = as.numeric(mcModel3[,"beta"]),
                                    dc = as.numeric(mcModel3[,"deltaC"]),
                                    ds = dsbiweekly,
                                    gamd = gamd)

  ciR0biweekly <- get95ci(cbind(mcR0biweekly))

  R0results <- c(R0wkly = meanR0weekly, R0wklylow = ciR0weekly[1], R0wklyhigh = ciR0weekly[2],
                 R0biwkly = meanR0biweekly, R0biwklylow = ciR0biweekly[1], R0biwklyhigh = ciR0biweekly[2])

  meanThrsh <- surveillancethreshold(bet = as.numeric(meanModel3["beta"]),
                                     dc = as.numeric(meanModel3["deltaC"]),
                                     gamd = gamd)

  mcThrsh <- surveillancethreshold(bet = as.numeric(mcModel3[,"beta"]),
                                   dc = as.numeric(mcModel3[,"deltaC"]),
                                   gamd = gamd)

  ciThrsh <- get95ci(cbind(mcThrsh))

  threshweeksresults <- c(threshweeks = 1/meanThrsh,
                          threshweekslow = 1/ciThrsh[2],
                          threshweekshigh = 1/ciThrsh[1]) / 7 * 0.954 * 0.85

  c(gamdtable,R0results,threshweeksresults)
}

gamdval <- gam * c(1, 2, 5, 10, 50, 100, Inf)

table4 <- matrix(0,length(gamdval),11)

for(i in seq_along(gamdval)){
  table4row <- getTable4row(gamdval[i])
  table4[i,] <- table4row
}
colnames(table4) <- names(table4row)
print(table4)

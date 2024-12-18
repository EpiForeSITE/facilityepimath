# This code produces the results in Figure 2 of the manuscript:
# "Transmission thresholds for the spread of infections in healthcare facilities"

rm(list=ls())

gam <- 1/387
eps <- 0.5
bet <- 0.051048866
deltaC <- 0.008447451
px <- 0.57993299
rx <- 0.02849861
rg <- 0.17911358
k <- 5.73518945

MGFmixedgamma <- function(x, prob, rate, shape, deriv=0)
  sum(exp(log(prob)+lgamma(shape+deriv)-lgamma(shape)-shape*log(1-x/rate)-deriv*log(rate-x)))

getMu <- function(px, rx, rg, k) px/rx + (1-px)*k/rg

getSigsq <- function(px, rx, rg, k)
  2*px/rx^2 + (1-px)*k*(k+1)/rg^2 - (px/rx + (1-px)*k/rg)^2

getR0intervention <- function(px, rx, rg, k){

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob = c(px, 1-px),
                                            rate = c(rx, rg),
                                            shape = c(1,k), deriv)

  facilityR0(S = 0, C = rbind(c(-deltaC-gam,0),c(deltaC,0)), A = rbind(1,0),
             transm = bet*c(1,1-eps), initS = 1, mgf = mgf)
}
getR0intervention <- Vectorize(getR0intervention)

rxRand <- rx*(0.6+runif(1000)*0.8)
rgRand <- rg*(0.6+runif(1000)*0.8)
kRand <- k*(0.6+runif(1000)*0.8)
pxRand <- px*(0.6+runif(1000)*0.8)

R0rand <- getR0intervention(rx=rxRand, rg=rgRand, k=kRand, px=pxRand)
muRand <- getMu(rx=rxRand, rg=rgRand, k=kRand, px=pxRand)
sigsqRand <- getSigsq(rx=rxRand, rg=rgRand, k=kRand, px=pxRand)

randPlot <- function(statRand, xlabel){
	plot(statRand, R0rand, pch='.', xlab=xlabel, ylab='Facility R0')
}

par(mfrow=c(2,2),mar=c(5,4,2,2)+.1)
randPlot(muRand, 'LOS mean')
randPlot(sqrt(sigsqRand), 'LOS standard deviation')
randPlot(sigsqRand/muRand, 'LOS variance to mean ratio (VMR)')
randPlot(muRand + sigsqRand/muRand, 'LOS mean plus VMR')

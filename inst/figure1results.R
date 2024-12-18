# This code produces the results in Figure 1 of the manuscript:
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

getMu <- function(px, rx, rg) px/rx + (1-px)*k/rg

getR0intervention <- function(px, rx, rg){

  mgf <- function(x, deriv=0) MGFmixedgamma(x, prob = c(px, 1-px),
                                            rate = c(rx, rg),
                                            shape = c(1,k), deriv)

  facilityR0(S = 0, C = rbind(c(-deltaC-gam,0),c(deltaC,0)), A = rbind(1,0),
             transm = bet*c(1,1-eps), initS = 1, mgf = mgf)
}
getR0intervention <- Vectorize(getR0intervention)

pts <- 100
chng <- seq(0,1,len=pts)
R0rxrg <- matrix(0,pts,1)
R0rx <- matrix(0,pts,1)
R0rg <- matrix(0,pts,1)
R0px <- matrix(0,pts,1)

murxrg <- matrix(0,pts,1)
murx <- matrix(0,pts,1)
murg <- matrix(0,pts,1)
mus <- matrix(0,pts,1)
mupx <- matrix(0,pts,1)

rxsBoth <- rx*(1+chng)
rgsBoth <- rg*(1+chng)
rxs <- rx*(1+chng)
rgs <- rg*(1+chng*1.5)
pxs <- px*(1-chng)

R0rxrg[,1] <- getR0intervention(px, rxsBoth, rgsBoth)
murxrg[,1] <- getMu(px, rxsBoth, rgsBoth)
R0rx[,1] <- getR0intervention(px, rxs, rg)
murx[,1] <- getMu(px, rxs, rg)
R0rg[,1] <- getR0intervention(px, rx, rgs)
murg[,1] <- getMu(px, rx, rgs)
R0px[,1] <- getR0intervention(pxs, rx, rg)
mupx[,1] <- getMu(pxs, rx, rg)

R0lim <- c(0.8,1.3); R0labs <- seq(0.8,1.3,0.1)
muLim <- c(25,35); muLabs <- seq(25,35,2)

getPlotVal <- function(val,lim) (val-lim[1])/diff(lim)

x <- matrix(rep(chng,1),pts,1)*100

par(mfrow=c(2,2))

doPlot <- function(R0vals,muVals,title){
  matplot(1/muVals,getPlotVal(R0vals,R0lim),lwd=2,
          type='l',lty=1,ylim=c(0,1),ylab = 'Facility R0',
          xlab = 'Mean length of stay', xlim = c(1/34, 1/26), axes=FALSE, main=title)
  box(); axis(1, 1/seq(34,26), seq(34,26))
  axis(2,getPlotVal(R0labs,R0lim),R0labs)
  lines(c(-100,200),getPlotVal(c(1,1),R0lim),col='grey',lty=3)
}

doPlot(R0rxrg,murxrg,'Increase discharge rate of all patients')
doPlot(R0rx,murx,'Increase discharge rate of high-variance patients')
doPlot(R0rg,murg,'Increase discharge rate of low-variance patients')
doPlot(R0px,mupx,'Decrease fraction of high-variance patients')

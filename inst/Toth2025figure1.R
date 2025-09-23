# This code produces the results in Figure 3 of the manuscript:
# "Transmission thresholds for the spread of infections in healthcare facilities" by Toth et al. (2025)

library(facilityepimath)

gam <- 1/387
eps <- 0.5
bet <- 0.051048866
dc <- 0.008447451
px <- 0.57993299
rx <- 0.02849861
rg <- 0.17911358
k <- 5.73518945

mgf <- function(x, deriv = 0) MGFmixedgamma(x, prob = c(px, 1 - px),
                                            rate = c(rx, rg),
                                            shape = c(1, k), deriv)

getInterventionEffect <- function(importRate, ds){
  S <- diag(0, 2)
  C <- rbind(c(-ds - dc - gam, 0, 0),
             c(ds, -dc - gam, 0),
             c(dc, dc, 0))
  A <- rbind(c(1, 0), c(0, 1 - eps), c(0, 0))
  transm <- bet * c(1 ,1 - eps, 1 - eps)
  initS <- rbind(1, 0)

  R0 <- facilityR0(S, C, A, transm, initS, mgf)

  R <- rbind(c(gam, 0, 0), c(0, gam, 0))
  init <- c((1 - importRate) * initS, importRate * c(1, 0, 0))

  eq <- facilityeq(S, C, A, R, transm, init, mgf)
  c(R0 = R0, clin = dc * sum(eq[3:4]))
}
getInterventionEffect <- Vectorize(getInterventionEffect)

snstvty <- 0.85
adherence <- 0.954
ds <- seq(0, 1 / 7, len = 100) * snstvty
dspts <- c(0, 1/(365 / 12), 1 / 14, 1 / 7) * snstvty * adherence
meanLOS <- mgf(0, 1)

out0001 <- getInterventionEffect(importRate = 0.0001, ds)
outpts0001 <- getInterventionEffect(importRate = 0.0001, dspts)

out001 <- getInterventionEffect(importRate = 0.001, ds)
outpts001 <- getInterventionEffect(importRate = 0.001, dspts)

out01 <- getInterventionEffect(importRate = 0.01, ds)
outpts01 <- getInterventionEffect(importRate = 0.01, dspts)

out1 <- getInterventionEffect(importRate = 0.05, ds)
outpts1 <- getInterventionEffect(importRate = 0.05, dspts)

oldpar <- par(mar = c(5, 4, 0.1, 2) + 0.1)
plot(out0001["R0", ], out0001["clin", ] * 1000 * meanLOS, xlab = expression(paste("Facility ", italic(R)[0])),
     ylab = 'Infections per 1000 admissions', type = 'l', lwd = 2, xlim = c(0.8, 1.3), ylim = c(0,71))
points(outpts0001["R0", ], outpts0001["clin", ] * 1000 * meanLOS, pch = 19)

lines(out001["R0", ], out001["clin", ] * 1000 * meanLOS, lwd = 2, lty=2)
points(outpts001["R0", ], outpts001["clin", ] * 1000 * meanLOS, pch = 19)

lines(out01["R0", ], out01["clin", ] * 1000 * meanLOS, lwd = 2, lty=3)
points(outpts01["R0", ], outpts01["clin", ] * 1000 * meanLOS, pch = 19)

lines(out1["R0", ], out1["clin", ] * 1000 * meanLOS, lwd = 2, lty=4)
points(outpts1["R0", ], outpts1["clin", ] * 1000 * meanLOS, pch = 19)

par(oldpar)

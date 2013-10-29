
library(foreach)
library(RColorBrewer)
cols <- brewer.pal(8,"Dark2")

n1 <- 15
n2 <- 15
delta <- 10.0
tau <- 10.0

no.sampled <- function(func, gamma1, gamma2) {
  samples <- times(100) %do% {
    func(n1, n2, delta, gamma1, gamma2, tau, 1000)
  }
  length(which(!is.na(samples)))
}


gammas <- seq(1e-9, 5, length.out = 100)
conf = expression(paste(n[1] == 15, ", ", n[2] == 15, ", ", 
                        delta == 10.0,  " and ", tau == 10.0))

opar <- par(mfrow=c(3,2))

## SAMPLER 1
y.0.1 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, g, 0.1))
y.1.0 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, g, 1.0))
y.1.9 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, g, 1.9))
y.2.5 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, g, 2.5))

plot(gammas, y.0.1, type='l', col=cols[1], ylim=c(0,100),
     main='First sampler',
     sub=conf,
     xlab = expression(gamma[1]), ylab = '% accepted')
lines(gammas, y.1.0, col=cols[2])
lines(gammas, y.1.9, col=cols[3])
lines(gammas, y.2.5, col=cols[4])

legend('bottomright', fill=cols[1:4],
       legend=c(expression(gamma[2] == 0.1),
                expression(gamma[2] == 1.0),
                expression(gamma[2] == 1.9),
                expression(gamma[2] == 2.5)))

y.0.1 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, 0.1, g))
y.1.0 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, 1.0, g))
y.1.8 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, 1.8, g))
y.1.9 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, 1.9, g))
y.2.0 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, 2.0, g))
y.2.5 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA, 2.5, g))

plot(gammas, y.0.1, type='l', col=cols[1], ylim=c(0,100),
     main='First sampler',
     sub=conf,
     xlab = expression(gamma[2]), ylab = '% accepted')
lines(gammas, y.1.0, col=cols[2])
lines(gammas, y.1.8, col=cols[3])
lines(gammas, y.1.9, col=cols[4])
lines(gammas, y.2.0, col=cols[5])
lines(gammas, y.2.5, col=cols[6])

legend('bottomright', fill=cols[1:6],
       legend=c(expression(gamma[1] == 0.1),
                expression(gamma[1] == 1.0),
                expression(gamma[1] == 1.8),
                expression(gamma[1] == 1.9),
                expression(gamma[1] == 2.0),
                expression(gamma[1] == 2.5)))



## SAMPLER 2
y.0.1 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, g, 0.1))
y.0.5 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, g, 0.5))
y.1.0 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, g, 1.0))
y.1.9 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, g, 1.9))
y.2.5 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, g, 2.5))

plot(gammas, y.0.1, type='l', col=cols[1], ylim=c(0,100),
     main='Second sampler',
     sub=conf,
     xlab = expression(gamma[1]), ylab = '% accepted')
lines(gammas, y.0.5, col=cols[2])
lines(gammas, y.1.0, col=cols[3])
lines(gammas, y.1.9, col=cols[4])
lines(gammas, y.2.5, col=cols[5])

legend('bottomright', fill=cols[1:5],
       legend=c(expression(gamma[2] == 0.1),
                expression(gamma[2] == 0.5),
                expression(gamma[2] == 1.0),
                expression(gamma[2] == 1.9),
                expression(gamma[2] == 2.5)))

y.0.1 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, 0.1, g))
y.1.0 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, 1.0, g))
y.1.9 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, 1.9, g))
y.2.5 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_alt, 2.5, g))

plot(gammas, y.0.1, type='l', col = cols[1], ylim=c(0,100),
     main='Second sampler',
     sub=conf,
     xlab = expression(gamma[2]), ylab = '% accepted')
lines(gammas, y.1.0, col=cols[2])
lines(gammas, y.1.9, col=cols[3])
lines(gammas, y.2.5, col=cols[4])

legend('bottomright', fill=cols[1:4],
       legend=c(expression(gamma[1] == 0.1),
                expression(gamma[1] == 1.0),
                expression(gamma[1] == 1.9),
                expression(gamma[1] == 2.5)))



## SAMPLER HYBRID
y.0.1 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, g, 0.1))
y.1.0 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, g, 1.0))
y.1.9 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, g, 1.9))
y.2.5 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, g, 2.5))

plot(gammas, y.0.1, type='l', col=cols[1], ylim=c(0,100),
     main='Hybrid sampler',
     sub=conf,
     xlab = expression(gamma[1]), ylab = '% accepted')
lines(gammas, y.1.0, col=cols[2])
lines(gammas, y.1.9, col=cols[3])
lines(gammas, y.2.5, col=cols[4])

legend('bottomright', fill=cols[1:4],
       legend=c(expression(gamma[2] == 0.1),
                expression(gamma[2] == 1.0),
                expression(gamma[2] == 1.9),
                expression(gamma[2] == 2.5)))

y.0.1 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, 0.1, g))
y.1.0 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, 1.0, g))
y.1.8 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, 1.8, g))
y.1.9 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, 1.9, g))
y.2.0 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, 2.0, g))
y.2.5 <- sapply(gammas, function(g) no.sampled(twoTumors_CaseB_TA_hybrid, 2.5, g))

plot(gammas, y.0.1, type='l', col=cols[1], ylim=c(0,100),
     main='Hybrid sampler',
     sub=conf,
     xlab = expression(gamma[2]), ylab = '% accepted')
lines(gammas, y.1.0, col=cols[2])
lines(gammas, y.1.8, col=cols[3])
lines(gammas, y.1.9, col=cols[4])
lines(gammas, y.2.0, col=cols[5])
lines(gammas, y.2.5, col=cols[6])

legend('bottomright', fill=cols[1:6],
       legend=c(expression(gamma[1] == 0.1),
                expression(gamma[1] == 1.0),
                expression(gamma[1] == 1.8),
                expression(gamma[1] == 1.9),
                expression(gamma[1] == 2.0),
                expression(gamma[1] == 2.5)))


par(opar)
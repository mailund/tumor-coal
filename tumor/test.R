
library(foreach)

n1 <- 10
n2 <- 15
delta <- 1.0
gamma1 <- 0.1
gamma2 <- 0.1
tau <- 1.0

TAs1 <- times(1000) %do% {
  twoTumors_CaseB_TA(n1, n2, delta, gamma1, gamma2, tau)
}

TAs2 <- times(1000) %do% {
  twoTumors_CaseB_TA_alt(n1, n2, delta, gamma1, gamma2, tau)
}

xlim <- range(TAs1, TAs2)
  
opar <- par(mfrow=c(2,1))
hist(TAs1, xlim=xlim)
hist(TAs2, xlim=xlim)
par(opar)

coaltimes.single.conditional <- function(n, delta, gamma, t1) {
  x1 <- gamma*exp(-delta*t1) / (1 - exp(-delta*t1))
  Us <- rev(sort(runif(n-1)))
  Zs <- Us / (1 + x1)
  Ts <- log( (1-Zs) / (1 - (1-gamma)*Zs)) / delta
  Ts
}

delta <- 1.0 ; t1 <- 1
gamma <- seq(0.1,0.9, length.out = 100)
y <- gamma*exp(-delta*t1) / (1 - exp(-delta*t1))
plot(gamma, y, type='l')

accept.prob.TA <- function(TA, n2, gamma2, delta) {
  if (gamma2 == 1) {
    x <- 1/n2
  } else {
    x <- -n2*gamma2 + sqrt((n2*gamma2)**2 + 4*(1-gamma2)) / (2*(1-gamma2))
  }
  
  M <- x*(1-x)**(n2-1) / (1-(1-gamma2)*x)**(n2+1)
  s <- delta * TA # TA = ttilde + tau
  
  (1/M*exp(-s)*(1-exp(-s))**(n2-1)) / ((1-(1-gamma2)*exp(-s))**(n2+1))
  
}

gamma <- seq(0.1, 1.9, length.out = 100)
probs <- sapply(gamma, function (g) accept.prob.TA(1, 15, g, 1))
plot(gamma, probs, type='l')
abline(h=0.0, lty='dashed')

gamma <- 0.5
Z <- seq(0,1,length.out = 100)
plot(Z, (1-Z)/(1-(1-gamma)*Z), type='l')
lines(Z, (1-Z)/(1-(1-0.1)*Z), col='red')
lines(Z, (1-Z)/(1-(1-0.9)*Z), col='blue')

n2 <- 5
gamma2 <- seq(0.5,1.5, length.out = 100)#seq(1.0001,3.0, length.out = 100)

opar <- par(mfrow=c(1,2))

y1 <- -n2 * gamma2
y2 <- sqrt((n2*gamma2)*(n2*gamma2) + 4*(1-gamma2)) / (2*(1-gamma2))

yli <- range(y1, y2)
plot(gamma2, y1, ylim=yli, type='l', col='blue', ylab='', xlab=expression(gamma[2]))
lines(gamma2, y2, col='red')
legend('bottomleft', fill=c('blue','red'),
       legend=c(expression(-n[2] * gamma[2]),
                expression(sqrt((n[2]*gamma[2])*(n[2]*gamma[2]) + 4*(1-gamma[2])) / (2*(1-gamma[2])))))

x <- -n2 * gamma2 + sqrt((n2*gamma2)*(n2*gamma2) + 4*(1-gamma2)) / (2*(1-gamma2))
plot(gamma2, x, type='l', xlab=expression(gamma[2]), ylab='x')
abline(h=0.0, lty='dashed')

M <- x*(1-x)**(n2-1) / (1-(1-gamma)*x)**(n2+1)

par(opar)
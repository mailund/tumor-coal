
library(ape)
opar <- par(mfrow=c(2,2))
tree <- rtumortree(15, 1, 1)
plot(tree)
plot(mutate(tree, 0.01))
plot(mutate(tree, 0.05))
plot(mutate(tree, 0.1))
par(opar)
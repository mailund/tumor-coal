
tsingle <- rtumourtree.single(4, 1, 1)
plot(mutate(tsingle,3))

tree <- rtumourtree.two(3, 3, 1, 1, 1, 1.2)
opar <- par(mfrow=c(2,1))
plot(tree)
plot(mutate(tree, 1))
par(opar)


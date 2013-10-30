
tau <- 1.2
tree <- rtumortree.two(5, 10, 1, 10, 10, tau)
plot(tree)
abline(v=c(0,tree$TA), col='red', lty='dashed')
abline(v=tree$TA - tau, col='blue', lty='dashed')

# If done correctly, these should all be TA
tree$TA
branching.times(tree$tree2)[1] + tree$tree2$root.edge
branching.times(tree$tree1)[1] + tree$tree1$root.edge + tau

#' Simulate a tumor tree for two tumours (case B).
#' 
#' Simulate a coalescence tree for case B in Wiuf's note. This is the genealogy
#' of two single tumours sharing a single driver mutation.
#' 
#' The time for the divergence of the two tumours is sampled as T_A in Wiuf's note, and
#' the two tumours are then sampled as single tumours in a similar way to 
#' \code{coaltimes.single}, but offset by parameter tau. See \code{coaltimes.single}
#' for details on how the individual tumour coalescence times are sampled.
#' 
#' @param n1      Number of samples / coalescence times in the first tumour.
#' @param n2      Number of samples / coalescence times in the second tumour.
#' @param delta   Parameter for the birth/death process.
#' @param gamma1  Parameter for the birth/death process in the first tumour.
#' @param gamma2  Parameter for the birth/death process in the second tumour.
#' @param tau     Offset in time of the sampling of the two tumours.
#' 
#' @export
rtumortree.two <- function(n1, n2, delta, gamma1, gamma2, tau) {
  # Get time points from the birth/death process in tumour 1
  t1s <- coaltimes.single(n1, delta, gamma1)
  # and then do the same for tumour 2
  t2s <- coaltimes.single(n2, delta, gamma2) + tau # FIXME: not sure tau should be added
  
  # Get the split time between the two tumours
  TA <- twoTumors_CaseB_TA(n1, n2, delta, gamma1, gamma2, tau)
  
  # Simulate two coalescence trees as single tumours, with the respective coal times.
  tree1 <- rcoal(n1, br = t1s[-1])
  tree2 <- rcoal(n2, br = t2s[-1])
  
  opar <- par(mfrow=c(2,1))
  plot(tree1, root.edge=TRUE)
  plot(tree2, root.edge=TRUE)
  par(opar)
  
  # Adjust the root edge time to TA
  tree1$root.edge <- TA - branching.times(tree1)[1]
  tree2$root.edge <- TA - branching.times(tree2)[1]
  
  opar <- par(mfrow=c(2,1))
  plot(tree1, root.edge=TRUE)
  plot(tree2, root.edge=TRUE)
  par(opar)
  
  # Merge the two trees
  tree <- tree1 + tree2
  
  tree$n1 <- n1
  tree$n2 <- n2
  tree$TA <- TA
  tree$t1s <- t1s
  tree$t2s <- t2s
  
  # Define as a sub-class of "phylo" and return the new object
  class(tree) <- c('tumortree.two.B', 'tumortree', 'phylo')
  tree
}


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
  # Get the split time between the two tumours
  TA <- twoTumors_CaseB_TA_hybrid(n1, n2, delta, 
                                  gamma1, gamma2, tau, 1000)
  if (is.na(TA))
    return (NA)

  # Get time points from the birth/death process in tumour 1
  t1s <- coaltimes.single.conditional(n1, delta, gamma1, TA)
  # and then do the same for tumour 2
  t2s <- coaltimes.single.conditional(n2, delta, gamma2, TA)
  
  # Simulate two coalescence trees as single tumours, with the respective coal times.
  tree1 <- rcoal(n1, br = coaltimes.differences(rev(t1s)))
  tree2 <- rcoal(n2, br = coaltimes.differences(rev(t2s)))
    
  # Adjust the root edge time to TA
  tree1$root.edge <- TA - branching.times(tree1)[1] - tau # FIXME?
  tree2$root.edge <- TA - branching.times(tree2)[1]
  
  # Merge the two trees
  tree <- tree1 + tree2
  
  tree$tree1 <- tree1
  tree$tree2 <- tree2
  tree$n1 <- n1
  tree$n2 <- n2
  tree$TA <- TA
  tree$t1s <- t1s
  tree$t2s <- t2s
  
  # Define as a sub-class of "phylo" and return the new object
  class(tree) <- c('tumortree.two.B', 'tumortree', 'phylo')
  tree
}

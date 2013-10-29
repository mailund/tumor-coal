
#' Simulate a tumor tree for a single tumor.
#' 
#' Simulate a coalescence tree for case A in Wiuf's note. This is the genealogy
#' of a single tumor where \code{n} cells were sampled.
#' 
#' Coalescence times are simulated using \code{coaltimes.single} and the parameters
#' of this function are the same as for that function. See the description of that function
#' for a description of the parameters. The tree topology is random.
#'
#' @param n      The number of coalescence times / leaves.
#' @param delta  Parameter for the birth/death process.
#' @param gamma  Parameter for the birth/death process.
#' 
#' @export
rtumortree.single <- function(n, delta, gamma) {
  # Get time points from the birht/death process
  ts <- coaltimes.single(n, delta, gamma)
  
  # Coalescence times are the birth times from t2 and down to tn
  coal.times <- ts[-1]
  
  # Simulate a coalescence tree with those times
  tree <- rcoal(n, br = coal.times)
  tree$root.edge <- ts[1] - ts[2] # t1 on the root edge, so the root edge has length t1 - t2
  tree$ts <- ts
  
  # Define as a sub-class of "phylo" and return the new object
  class(tree) <- c('tumortree.single', 'tumortree', 'phylo')
  tree
}


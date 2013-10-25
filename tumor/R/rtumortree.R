birth.times <- function(n, delta, gamma) {
  Us <- sort(runif(n))
  -1.0/delta * log( Us / (gamma + (1-gamma) * Us) )
}

#' Simulate a tumor tree.
#' 
#' @export
rtumortree <- function(n, delta, gamma) {
  # Get time points from the birht/death process
  ts <- birth.times(n, delta, gamma)
  
  # Coalescence times are the birth times from t2 and down to tn
  coal.times <- ts[-1]
  
  # Simulate a coalescence tree with those times
  tree <- rcoal(n, br = coal.times)
  tree$ts <- ts
  
  # Define as a sub-class of "phylo" and return the new object
  class(tree) <- c('tumortree', 'phylo')
  tree
}


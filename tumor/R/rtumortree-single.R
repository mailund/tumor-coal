
#' Simulate the coalescence times for a single tumor.
#' 
#' Simulate coalescence times for a tumor based on a birth/death process
#' with birth rate lambda and death rate mu. The simulation uses
#' scaled parameters delta = lambda - mu and gamma = (lambda-mu)/(lambda*rho)
#' where rho is the fraction of tumor cells sequenced.
#' 
#' See Wiuf's note for details.
#' 
#' The function will return \code{n} time points, where t1 is the time for the 
#' driver mutation and the remaining are the actual coalescence times for samples
#' from the tumor.
#' 
#' @param n      The number of coalescence times.
#' @param delta  Parameter for the birth/death process.
#' @param gamma  Parameter for the birth/death process.
#' 
#' @export
coaltimes.single <- function(n, delta, gamma) {
  Us <- sort(runif(n))
  -1.0/delta * log( Us / (gamma + (1-gamma) * Us) )
}

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
  tree$ts <- ts
  
  # Define as a sub-class of "phylo" and return the new object
  class(tree) <- c('tumortree.single', 'tumortree', 'phylo')
  tree
}


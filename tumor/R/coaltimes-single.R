
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

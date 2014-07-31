
sample.times.B1 <- function(n1, n2, delta, gamma1, gamma2, tau) {
  
  # Functions for rejection sampling
  R <- function(z, t2, s2, gamma1.tilde, gamma2.tilde) {
    nom <- (delta * gamma2.tilde * exp( -delta*(t2-s2) ) * z * (gamma1.tilde - (gamma1.tilde - 1)*z))
    denom <- ( gamma1.tilde + (gamma2.tilde - 1)*exp(-delta*(t2-s2) - (gamma1.tilde-1))*z )**2
    nom / denom
  }
  gamma1.tilde <- function(t2) { 1 + (gamma1 - 1) * exp( -delta*(t2-tau) ) }
  gamma2.tilde <- function(s2) { 1 + (gamma2 - 1) * exp( -delta*s2 ) }
  
  R.max <- function(t2, s2, gamma1.tilde, gamma2.tilde) {
    in.unit.interval <- function(x) x > 0 && x < 1
    z.tilde <- gamma1.tilde / (gamma1.tilde - 1 + (gamma2.tilde-1)*exp(-delta*(t2-s2)))
    z.max <- ifelse(in.unit.interval(z.tilde), z.tilde, 1)
    ifelse(in.unit.interval(z.max),
           delta*gamma2.tilde / (4*(gamma2.tilde - 1)),
           delta*gamma2.tilde*exp(-delta*(t2-s2)) / (1+(gamma2.tilde-1)*exp(-delta*(t2-s2)))**2
    )
  }
  
  acceptance.prob <- function(z, t2, s2, gamma1.tilde, gamma2.tilde) {
    R(z, t2, s2, gamma1.tilde, gamma2.tilde) / R.max(t2, s2, gamma1.tilde, gamma2.tilde)
  }
  
  # Sample times and reject/accept...
  while(TRUE) {
    t1s <- coaltimes.single(n1, delta, gamma1) + tau
    t2s <- coaltimes.single(n2, delta, gamma2)
  
    t2 <- t1s[2]
    s2 <- t2s[2]
    g1.tilde <- gamma1.tilde(t2)
    g2.tilde <- gamma2.tilde(s2)
    
    z <- runif(1)
    
    if (t2 > s2) {
      A <- acceptance.prob(z, t2, s2, g1.tilde, g2.tilde)
      p <- runif(1)
      cat("A =", A, "p =", p, "\n")
      if (p <= A) {
        TA <- -log(z/(g1.tilde - (g1.tilde-1)*z))/delta + t2
        return(list(TA=TA, t1s=t1s, t2s=t2s))
      }
    } else {
      A <- acceptance.prob(z, s2, t2, g2.tilde, g1.tilde)
      p <- runif(1)
      cat("A =", A, "p =", p, "\n")
      if (p <= A) {
        TA <- -log(z/(g2.tilde - (g2.tilde-1)*z))/delta + s2
        return(list(TA=TA, t1s=t1s, t2s=t2s))
      }
    }
  }
}

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
rtumourtree.two <- function(n1, n2, delta, gamma1, gamma2, tau) {

  times <- sample.times.B1(n1, n2, delta, gamma1, gamma2, tau)
  t1s <- times$t1s
  t2s <- times$t2s
  TA <- times$TA
  
  # Simulate two coalescence trees as single tumours, with the respective coal times.
  tree1 <- rcoal(n1, br = coaltimes.differences(rev(t1s - tau)))
  tree2 <- rcoal(n2, br = coaltimes.differences(rev(t2s)))
  
  # Adjust the root edge time to TA
  tree1$root.edge <- TA - branching.times(tree1)[1] - tau
  tree2$root.edge <- TA - branching.times(tree2)[1]
  
  # Merge the two trees
  tree <- tree1 + tree2
  
  tree$tree1 <- tree1
  tree$tree2 <- tree2
  tree$n1 <- n1
  tree$n2 <- n2
  tree$TA <- TA
  tree$tau <- tau
  tree$t1s <- t1s
  tree$t2s <- t2s
  
  # Define as a sub-class of "phylo" and return the new object
  class(tree) <- c('tumourtree.two.B', 'tumourtree', 'phylo')
  tree
}


library(ape)

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

#' Place mutations on a tumor tree.
#' 
#' @export
mutate <- function(tree, theta) {
  # Place mutations on the tree based on mutation rate theta and edge lengths.
  lambdas <- tree$edge.length * theta
  no.muts <- rpois(length(tree$edge.length), lambdas)
  tree$no.muts <- no.muts
  tree
}

# Helper function for adding mutations to a plot of a tree
plot.tumortree.mutations <- function(tree) {
  for (edge in 1:length(tree$no.muts)) {
    k <- tree$no.muts[edge]
    if (k > 0) {
      l <- tree$edge.length[edge]
      adjs <- runif(k, min=-l/5.0, max=l/2.0) # range based on eyeballed aestetics
      for (i in 1:k) {
        edgelabels('', frame='none', col='black', bg='red', pch=21, edge=edge, adj=adjs[i])
      }
    }
  }
}

#' Plot a tumor tree
#' 
#' Plots the tree (using \code{ape} plotting) including mutations if the tree includes those.
#' 
#' @param tree The tree to plot
#' @param ...  Options passed to \code{plot.phylo}
#' 
#' @export
plot.tumortree <- function(tree, ...) {
  plot.phylo(tree, ...)
  if (!is.null(tree$no.muts)) {
    plot.tumortree.mutations(tree)
  }
  axisPhylo()
}




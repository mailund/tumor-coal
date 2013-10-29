
# This file contains polymorphic code for the class "tumortree".


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




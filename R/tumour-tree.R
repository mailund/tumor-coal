
# This file contains polymorphic code for the class "tumourtree".


# Helper function for adding mutations to a plot of a tree
plot.tumourtree.mutations <- function(tree) {
  for (edge in 1:length(tree$no.muts)) {
    k <- tree$no.muts[edge]
    if (k > 0) {
      l <- tree$edge.length[edge]
      adjs <- runif(k, min=1/2 - l/2, max=1/2 + l/2)
      for (i in 1:k) {
        edgelabels('', frame='none', col='black', bg='red', pch=21, edge=edge, adj=adjs[i])
      }
    }
  }
}

#' Plot a tumour tree
#' 
#' Plots the tree (using \code{ape} plotting) including mutations if the tree includes those.
#' 
#' @param tree The tree to plot
#' @param ...  Options passed to \code{plot.phylo}
#' 
#' @export
plot.tumourtree <- function(tree, ...) {
  plt <- plot.phylo(tree, ...)
  if (!is.null(tree$no.muts)) {
    plot.tumourtree.mutations(tree)
  }
  axisPhylo()
  invisible(plt)
}

#' Plot a tumor tree
#' 
#' Plots the tree (using \code{ape} plotting) including mutations if the tree includes those.
#' 
#' @param tree The tree to plot
#' @param ...  Options passed to \code{plot.phylo}
#' 
#' @export
plot.tumourtree.single <- function(tree, ...) {
  plt <- plot.tumourtree(tree, root.edge=TRUE, ...)
  abline(v=0, col='black', lty='dashed')
  
  n <- length(tree$ts)
  scaled_times <- -(tree$ts - tree$ts[1]) # The ape package does weird stuff with the x-axis
  
  abline(v = scaled_times[1], col='red', lty='dashed')
  mtext(substitute(t[i], list(i=1)), at = scaled_times[1], col='red')
  mtext("Cancer\norigin", line=2, at = scaled_times[1], col='red')
  
  abline(v = scaled_times[2], col='black', lty='dashed')
  mtext(substitute(t[i], list(i=2)), at = scaled_times[2], col='black')
  mtext("MRCA", line=2, at = scaled_times[2], col='black')
  
  for (i in 3:n) {
    abline(v = scaled_times[i], col='blue', lty='dashed')
    mtext(substitute(t[i], list(i=i)), at = scaled_times[i], col='blue')  
  }
  
  invisible(plt)
}

#' Plot a tumor tree
#' 
#' Plots the tree (using \code{ape} plotting) including mutations if the tree includes those.
#' 
#' @param tree The tree to plot
#' @param ...  Options passed to \code{plot.phylo}
#' 
#' @export
plot.tumourtree.two.B <- function(tree, ...) {
  plt <- plot.tumourtree(tree, ...)
  abline(v=0, col='blue', lty='dashed')
  mtext(expression(t[A]), at=0, col='blue')
  abline(v=tree$TA - tree$tau, col='blue', lty='dashed')
  mtext(expression(tau), at=tree$TA - tree$tau, col='blue')
  invisible(plt)
}



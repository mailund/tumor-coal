
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
  plt <- plot.phylo(tree, show.tip.label=FALSE, ...)
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
  
  rightmost <- max(node.depth.edgelength(tree))
  rescale <- function(x) -(x - rightmost)
  
  plt <- plot.tumourtree(tree, ...)
  
  abline(v=0, col='black', lty='dashed')
  mtext(expression(t[A]), line=1, at=0, col='black')
  
  abline(v=rescale(tree$tau), col='black', lty='dashed')
  mtext(expression(tau), line=1, at=rescale(tree$tau), col='black')
  abline(v=rightmost, col='black', lty='dashed')
  mtext(0, line=1, at=rightmost, col='black')
  
  for (i in 2:tree$n1) {
    abline(v = rescale(tree$t1s[i]), col='blue', lty='dotted')
    mtext(substitute(t[i], list(i=i)), at = rescale(tree$t1s[i]), col='blue')  
  }
  for (i in 2:tree$n2) {
    abline(v = rescale(tree$t2s[i]), col='darkgreen', lty='dotted')
    mtext(substitute(s[i], list(i=i)), at = rescale(tree$t2s[i]), col='darkgreen')  
  }
  
  invisible(plt)
}



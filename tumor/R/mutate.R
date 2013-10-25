
#' Place mutations on a tumor tree.
#' 
#' Mutations are placed on the branches of the tree by a Poison process, so each branch of length \code{L}
#' gets a number of mutations drawn as \code{prois(L*theta)} where \code{theta} is the (scaled)
#' mutation rate.
#' 
#' @param tree   The tree to add mutations to.
#' @param theta  The scaled mutation rate.
#' 
#' @export
mutate <- function(tree, theta) {
  # Place mutations on the tree based on mutation rate theta and edge lengths.
  lambdas <- tree$edge.length * theta
  no.muts <- rpois(length(tree$edge.length), lambdas)
  tree$no.muts <- no.muts
  tree
}
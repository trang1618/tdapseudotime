#' Calculate the pairwise similarity between rows of numeric matrices
#'
#' @param a Input matrix.
#' @param b Input matrix. Optional.
#' If `b` is NULL, use `a`, i.e. similarity between the rows of `a`
#' would be computed.
#'
#' @return Matrix of the similarity between the rows of a and
#' the rows of b.
#'
#' @export
#'
#' @examples
#' m <- matrix(1:8, ncol = 4) # a 2x4 matrix
#' n <- matrix(9:16, ncol = 4) # a 2x4 matrix
#' cosine_sim_func(m)
#' cosine_sim_func(m, n) # a 4x4 similarity matrix
cosine_sim_func <- function(a, b = NULL){
  a_norm <- a / sqrt(rowSums(a * a))
  if (is.null(b)){
    b_norm <- a_norm
  } else {
    b_norm <- b / sqrt(rowSums(b * b))
  }
  cosine_sim <- a_norm %*% t(b_norm)
  cosine_sim[cosine_sim > 1] <- 1.0
  cosine_sim
}

# ==============================================================
# SIMILARITY FUNCTIONS

#' Compute Jaccard similarity to assign each subject to
#' most similar trajectory
#'
#' @param x Numeric vector of node numbers
#' indicating an individual's trajectory
#' @param y Numeric vector of node numbers
#' indicating a general trajectory
#'
#' @return A scalar as measure of similarity
#' @export
#'
#' @examples
#' sim_jaccard(c(1,2), c(1,3,4))
sim_jaccard  <- function(x, y) {
  (length(intersect(x, y))) / length(union(x, y))
}

#' Compute intersection similarity to assign each subject to
#' most similar trajectory
#'
#' @param x Numeric vector of node numbers
#' indicating an individual's trajectory
#' @param y Numeric vector of node numbers
#' indicating a general trajectory
#'
#' @return A scalar as measure of similarity
#' @export
#'
#' @examples
#' sim_intersection(c(1,2), c(1,3,4))
sim_intersection  <- function(x, y) {
  length(intersect(x, y)) / length(x)
}

#' Compute similarity based on length to assign each subject to
#' most similar trajectory
#'
#' @param x Numeric vector of node numbers
#' indicating an individual's trajectory
#' @param y Numeric vector of node numbers
#' indicating a general trajectory
#'
#' @return A scalar as measure of similarity
#' @export
#'
#' @examples
#' sim_length(c(1,2), c(1,3,4))
sim_length  <- function(x, y) {
  exp(-abs(length(x) - length(y)))
}


#' Find the shortest paths.
#'
#' @param x Pair of node ids from start node to end node.
#' @param mst_weights Output from igraph::mst.
#'
#' @return Vector of node ids of shortest paths.
#'
shortest_paths_func <- function(x, mst_weights){
  igraph::all_shortest_paths(
    mst_weights,
    from = x[1],
    to = x[2],
    mode = c("out"),
    weights = NULL
  )
}

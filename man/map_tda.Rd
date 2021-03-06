% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-funcs.R
\name{map_tda}
\alias{map_tda}
\title{Run TDAMapper on a matrix of lab values}
\usage{
map_tda(lab_values_mat, ...)
}
\arguments{
\item{lab_values_mat}{Matrix of (imputed) lab values.}

\item{\dots}{Arguments to pass to TDAMapper()}
}
\value{
TDAMapper object which is a list of items named \code{adjacency}
(adjacency matrix for the edges), \code{num_vertices} (integer number of vertices),
\code{level_of_vertex} (vector with \code{level_of_vertex[i]} = index of
the level set for vertex i), \code{points_in_vertex} (list with
\code{points_in_vertex[[i]]} = vector of indices of points in vertex i),
\code{points_in_level} (list with \code{points_in_level[[i]]} = vector of indices
of points in level set i, and \code{vertices_in_level} (list with
\code{vertices_in_level[[i]]} = vector of indices of vertices in level set i.
}
\description{
Run TDAMapper on a matrix of lab values
}
\examples{
my_tda <- map_tda(scaled_lab_mat)
str(my_tda, max.lev = 1)
}

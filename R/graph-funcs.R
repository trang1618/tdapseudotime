#' Run TDAMapper on a matrix of lab values
#'
#' @param lab_values_mat Matrix of (imputed) lab values.
#' @param \dots Arguments to pass to TDAMapper()
#'
#' @return TDAMapper object which is a list of items named \code{adjacency}
#' (adjacency matrix for the edges), \code{num_vertices} (integer number of vertices),
#' \code{level_of_vertex} (vector with \code{level_of_vertex[i]} = index of
#' the level set for vertex i), \code{points_in_vertex} (list with
#' \code{points_in_vertex[[i]]} = vector of indices of points in vertex i),
#' \code{points_in_level} (list with \code{points_in_level[[i]]} = vector of indices
#' of points in level set i, and \code{vertices_in_level} (list with
#' \code{vertices_in_level[[i]]} = vector of indices of vertices in level set i.
#'
#' @export
#'
#' @examples
#' my_tda <- map_tda(scaled_lab_mat)
#' str(my_tda, max.lev = 1)
map_tda <- function(lab_values_mat, ...){
  # Preparing for running TDA Mapper
  # Compute cosine similarity and principal and secondary SVD
  cosine_sim <- cosine_sim_func(lab_values_mat)
  trunc_svds <- RSpectra::svds(cosine_sim, k = 2, nu = 2, nv = 2)
  svd1 <- - trunc_svds$u[, 1]
  svd2 <- - trunc_svds$u[, 2]

  f_sim_map <- TDAmapper::mapper2D(
    distance_matrix = cosine_sim,
    filter_values = list(svd1 , svd2),
    ...
  )

  f_sim_map
}

#' Make igraph object by calculating adjacency, resizing nodes,
#' weighting edges and color nodes.
#'
#' @param f_sim_map TDAmapper object
#' @param data Processe data frame of the original data.
#' @param enrich_var Character string of enrichment variable.
#' Any column name of the processed data frame.
#' @param color_method Character string specifying the coloring method.
#' Can be 'basic', 'clust_shade' or 'clust_color'.
#' @param my_colors Character vector of hex values specifying
#' color palette for enrichment.
#'
#' @importFrom dplyr left_join arrange mutate filter bind_rows distinct select pull
#' @return igraph object of the graph output.
#' @export
#'
#' @examples
#' my_tda <- map_tda(scaled_lab_mat)
#' make_tda_graph(my_tda, sim_dat, 'time')
#'
make_tda_graph <- function(
  f_sim_map, data, enrich_var, color_method = 'clust_color',
  my_colors = c("#00A3DD", "#60C659", "#FFBC21", "#FF7F1E", "#EF2B2D")) {
  f_time <- data[, c('id', enrich_var)] %>% `colnames<-`(c('ID', 'val'))
  f_graph <- igraph::graph.adjacency(f_sim_map$adjacency, mode = "undirected")
  f_graph <- resize_nodes(f_sim_map, f_graph)
  f_graph <- weight_edges(f_sim_map, f_graph, f_time)
  f_graph <- color_graph(f_sim_map, f_graph, f_time, my_colors, color_method)
  f_graph
}

#' Get minimum spanning tree from graph
#'
#' @param f_graph igraph object, out put from graph.adjacency
#'
#' @return igraph object representing the minimumspanningtree
#'
minspantree <- function(f_graph){
  tree <- igraph::mst(f_graph, weights = f_graph$clusters$edge.betweenness)
  plot(tree)
  graphics::legend(
    "topright",
    legend = f_graph$pal$cluster,
    col = f_graph$pal$color,
    fill = f_graph$pal$color,
    horiz = TRUE,
    cex = 0.4
  )
  tree
}


#' Get basic properties of graph
#'
#' @param f_graph igraph object, out put from graph.adjacency
#' @param f_sim_map Optional TDAmapper object
#' @param simplified Boolean. If TRUE, return a numeric vectors of
#' basic properties. If FALSE, return a more comprehensive list of
#' properties including maximal connected components and
#' number of elements in each node.
#'
#' @return List or numeric vector of different properties including
#' number of graph components, number of nodes, median degree,
#' edge density, etc.
#'
#' @export
#'
#' @examples
#' my_tda <- map_tda(scaled_lab_mat)
#' my_graph <- make_tda_graph(my_tda, sim_dat, 'time')
#' get_graph_properties(my_graph)
get_graph_properties <- function(f_graph, f_sim_map = NULL, simplified = TRUE){
  comps <- igraph::components(f_graph) # maximal connected components of a graph
  out_list <- list(
    n_comps = length(unique(comps$membership)),
    # number of nodes
    n_nodes = length(V(f_graph)),
    # number of observation in nodes
    median_degree = stats::median(igraph::degree(f_graph)),
    edge_density = igraph::edge_density(f_graph),
    clique_length = length(igraph::cliques(f_graph)),
    diameter = igraph::diameter(f_graph),
    mean_distance = igraph::mean_distance(f_graph)
  )
  if (simplified){
    unlist(out_list)
  } else {
    if (!is.null(f_sim_map)){
      n_elements_in_nodes <- sapply(f_sim_map$points_in_vertex, length)
    } else {
      n_elements_in_nodes <- NULL
    }
    c(out_list,
      list(comps = comps, n_elements_in_nodes = n_elements_in_nodes))
  }
}


#' Color nodes of graph given coloring method.
#'
#' @param f_sim_map TDAmapper object
#' @param f_graph igraph object, out put from graph.adjacency
#' @param f_time Data frame of the original data but
#' with only two columns: ID and val (e.g., time)
#' @param method Character string specifying the coloring method.
#' Can be 'basic', 'clust_shade' or 'clust_color'.
#' @param my_colors Character vector of hex values specifying
#' color palette for enrichment.
#'
#' @importFrom igraph V
#' @return Modified graph with colors at nodes.
#'
color_graph <- function(
  f_sim_map, f_graph, f_time, my_colors,
  method = c('basic', 'clust_shade', 'clust_color', 'none')) {

  method <- match.arg(method)

  if (method == 'none'){
    plot(f_graph)
  } else if (method == 'basic'){
    color_map <- color_vertex(f_sim_map, f_time, my_colors)
    igraph::V(f_graph)$color <- color_map$colors
    plot(f_graph)
  } else if (method == 'clust_shade'){
    my_clusters <- commu_clus(f_graph) # community detection
    plot(my_clusters, f_graph) # highlight community clusters
    # enrich graph with cluster colors
    igraph::V(f_graph)$color <- paste(my_clusters$membership)
    f_graph$clusters <- my_clusters
  } else {
    my_clusters <- commu_clus(f_graph) # community detection
    node_color <- color_clust(f_sim_map, my_clusters)
    igraph::V(f_graph)$color <- node_color$color
    pal <- node_color %>% select(- node) %>% distinct()
    plot(f_graph) # plot with assigned palette
    graphics::legend(
      "topleft",
      legend = pal$cluster,
      col = pal$color,
      fill = pal$color,
      horiz = TRUE,
      box.lty = 0,
      cex = 0.8
    )
    f_graph$node_color <- node_color
    f_graph$clusters <- my_clusters
    f_graph$pal <- pal
  }
  f_graph
}


#' Resize the nodes of graph given the number of points it contains
#'
#' @param f_sim_map TDAmapper object
#' @param f_graph igraph object, out put from graph.adjacency
#' @param mult Scalar as a multiplicative size rescale parameter.
#' @param add Scalar as an additive size rescale parameter.
#'
#' @return An igraph object of the modified f_graph.
#'
resize_nodes <- function(f_sim_map, f_graph, mult = 6, add = 8){
  n_vertices <- f_sim_map$num_vertices
  vertex.size <- vector('numeric', n_vertices)
  for (i in seq.int(n_vertices)) {
    points.in.vertex <- f_sim_map$points_in_vertex[[i]]
    vertex.size[i] <- length(f_sim_map$points_in_vertex[[i]])
  }
  min_size <- min(vertex.size)
  max_size <- max(vertex.size)
  igraph::V(f_graph)$size <- (vertex.size - min_size) / (max_size - min_size) * mult + add
  f_graph
}

#' Weight edges of a graph based on mean time of each edge.
#'
#' @param f_sim_map TDAmapper object
#' @param f_graph igraph object, out put from graph.adjacency
#' @param f_time Data frame of the original data but
#' with only two columns: ID and val (time)
#'
#' @return An igraph object of the modified f_graph.
#'
weight_edges <- function(f_sim_map, f_graph, f_time){
  n_edges <- length(igraph::E(f_graph))

  for (j in seq.int(n_edges)) {
    tail <- igraph::tail_of(f_graph, igraph::E(f_graph)[j])
    head <- igraph::head_of(f_graph, igraph::E(f_graph)[j])
    pointInTail <- f_sim_map$points_in_vertex[[tail]]
    pointInHead <- f_sim_map$points_in_vertex[[head]]
    commonIDS <- intersect(pointInTail, pointInHead)

    # ========  TIME WEIGHT
    igraph::E(f_graph)$weight[j] <- f_time %>%
      filter(ID %in% f_time$ID[commonIDS]) %>%
      pull(val) %>%
      mean()
  }
  f_graph
}

#' Adjust colors of nodes using enrichment function.
#'
#' @param f_sim_map TDAmapper object
#' @param f_time Data frame of the original data but
#' with only two columns: ID and val (time)
#' @param my_colors Character vector of hex values specifying
#' color palette for enrichment.
#'
#' @return Data frame of vertex ids and colors.
#'
color_vertex <- function(f_sim_map, f_time, my_colors){
  colfunc <- grDevices::colorRampPalette(my_colors)
  y.mean.vertex <- list()
  for (i in 1:f_sim_map$num_vertices) {
    points.in.vertex <- f_sim_map$points_in_vertex[[i]]
    y.mean.vertex[[i]] <-
      data.frame(id = paste(i),
                 value = f_time %>%
                   filter(ID %in% f_time$ID[points.in.vertex]) %>%
                   pull(val) %>%
                   as.numeric() %>%
                   mean())
  }

  color_map <- bind_rows(y.mean.vertex) %>%
    arrange(value) %>%
    mutate(colors = unique(.$value) %>% length() %>% colfunc()) %>%
    arrange(as.numeric(id))
    # check: order of the igraph::V(f_graph) vertices is the same of the clrMap$id
}

#' Detect community structure based on edge igraph::betweeness.
#'
#' @param f_graph igraph object, out put from graph.adjacency
#' @param directed Logical. whether to calculate directed edge betweenness
#' for directed graphs. Ignored for undirected graphs.
#' @param bridges Logical. whether to return a list the edge removals
#' which actually a component of the graph.
#' @param \dots Additional arguments to pass to igraph::edge.betweenness.community.
#'
#' @return Community clusters from igraph.
#'
commu_clus <- function(f_graph, directed = FALSE, bridges = TRUE, ...){
  igraph::edge.betweenness.community(
    f_graph,
    weights = igraph::E(f_graph)$value,
    directed = directed,
    bridges = bridges,
    ...
  )
}

#' Get cluster colors for nodes.
#'
#' @param f_sim_map TDAmapper object
#' @param my_clusters Community clusters from igraph.
#'
#' @return Data frame of nodes and corresponding colors
#' based on the cluster each node belongs.
#'
color_clust <- function(f_sim_map, my_clusters) {
  cluster_vec <- as.factor(unique(my_clusters$membership))
  # Make a palette of cluster colors
  my_palette <- data.frame(
    color = RColorBrewer::brewer.pal(length(cluster_vec), "Set1"),
    cluster = cluster_vec
  )

  # Create data frame of nodes and cluster
  node_color <- data.frame(
    node = f_sim_map$level_of_vertex,
    cluster = as.factor(my_clusters$membership)
  ) %>%
    left_join(my_palette, by = 'cluster') %>%
    arrange(node)
}



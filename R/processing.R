#' Assign observations (ids) to nodes in network.
#'
#' @param f_sim_map TDAmapper object
#' @param processed_data Dataframe. Processed from original data.
#'
#' @import TDAmapper
#' @importFrom igraph V E mst edge.betweenness.community betweenness graph.adjacency
#' @return Dataframe of id, node and covid_id.
#'
convert_id_to_node <- function(f_sim_map, processed_data) {
  rowid_to_nodes <- f_sim_map$points_in_vertex
  # print(names(rowid_to_nodes))
  if (is.null(names(rowid_to_nodes)))
    names(rowid_to_nodes) <- seq.int(length(rowid_to_nodes))
  node_reps <- rep(names(rowid_to_nodes), lengths(rowid_to_nodes))
  # # Check if this this is the desired vector
  # table(node_reps)
  # lengths(rowid_to_nodes)

  rowid_node_df <- data.frame(
    id = unlist(rowid_to_nodes),
    node = as.integer(node_reps))
}

#' Find trajectories in the MST
#' All the trajectory or choose the nodes.
#'
#' @param processed_data Dataframe. Processed from original data.
#' @param f_sim_map TDAmapper object.
#' @param f_graph igraph object, output from graph.adjacency.
#'
#' @return List of trajectories computed from shortest paths.
#'
#' @export
#'
#' @examples
#' my_tda <- map_tda(scaled_lab_mat)
#' my_graph <- make_tda_graph(my_tda, sim_dat, 'time')
#' my_trajs <- find_trajectories(sim_dat, my_tda, my_graph)
#' head(my_trajs[[1]])
#' head(my_trajs[[2]])
#'
find_trajectories <- function(processed_data, f_sim_map, f_graph) {
  minspantreeweights <- minspantree(f_graph)
  mst_between <- betweenness(minspantreeweights)
  starting_nodes <- V(minspantreeweights)[mst_between == min(mst_between)]
  ending_nodes <- V(minspantreeweights)[mst_between == max(mst_between)]
  starts_and_ends <- expand.grid(starting_nodes, ending_nodes)
  #=========== OR chose specific nodes in the MST
  #starting_nodes<-c(3,4,8,20,23)
  #ending_nodes<-c(25,18)

  trajectories <- apply(
    starts_and_ends, 1, shortest_paths_func, mst_weights = minspantreeweights)
  trajectories
}

#' Compute similarity based on trajectory
#'
#' @param processed_data Dataframe. Processed from original data.
#' @param f_graph igraph object, out put from graph.adjacency.
#' @param trajectories List of trajectories computed from shortest paths.
#' @param f_sim_map TDAmapper object
#' @param verbose Boolean. Whether to print trajectory.
#'
#' @return List of two components:
#' dataframe with various similarity measures and
#' dataframe of the observations mapped to corresponding node and cluster.
#'
#' @export
#' @examples
#' my_tda <- map_tda(scaled_lab_mat)
#' my_graph <- make_tda_graph(my_tda, sim_dat, 'time')
#' my_trajs <- find_trajectories(sim_dat, my_tda, my_graph)
#' compute_similarity(sim_dat, my_graph, my_trajs, my_tda)
#'
compute_similarity <- function(
  processed_data, f_graph, trajectories, f_sim_map = NULL, verbose = TRUE) {

  if (!'node' %in% colnames(processed_data)){
    rowid_node_df <- convert_id_to_node(f_sim_map, processed_data)
    processed_data <- processed_data %>%
      left_join(rowid_node_df, by = 'id')
  }

  out_data <- processed_data %>%
    left_join(f_graph$node_color, by = 'node') %>%
    select(covid_id, id, node, cluster) %>%
    distinct()

  unique_patients <- unique(out_data$covid_id)
  similarity_ls <- list()

  for (i in unique_patients) {
    patienti_dat <- out_data %>% filter(covid_id == i)
    traj <- patienti_dat$node

    # Define trajectoris among clusters (coarse granularity)
    trajcluster <- unique(patienti_dat$cluster)

    temp <- list()

    for (t in 1:length(trajectories)) {
      if (length(trajectories[[t]]$res) > 0) {
        traj_res_1 <- trajectories[[t]]$res[[1]]
        temp[[t]] <- data.frame(
          covid_id = i,
          trajPaz = paste(traj, collapse = " "),
          trajPazclusters = paste(stats::na.omit(trajcluster), collapse = " "),
          trajNumb = t,
          trajElmnts = paste(traj_res_1, collapse = " "),
          trajLenght = length(traj_res_1),
          SJ = sim_jaccard(traj, traj_res_1),
          SI = sim_intersection(traj, traj_res_1),
          SL = sim_length(traj, traj_res_1),
          # Jaro-Winkler distance > take into account order
          # Values are comparable to Jaccard but formally is more correct (TBD)
          JW = stringdist::stringsim(
            paste(traj_res_1, collapse = "_"),
            paste(traj, collapse = "_"),
            method = 'jw'
          )
        )
      }
    }
    similarity_ls[[i]] <- bind_rows(temp)
  }
  traj_clusters <- add_clust_info(f_graph, trajectories)
  if (verbose) print(traj_clusters)
  similarity_df <- bind_rows(similarity_ls) %>%
    left_join(traj_clusters, by = 'trajElmnts')

  list(similarity_df, out_data)
}



#' Add cluster information
#'
#' @param f_graph igraph object, out put from graph.adjacency
#' @param trajectories List of trajectories computed from shortest paths.
#'
#' @return Dataframe of trajectory elements (trajElmnts)
#' and cluster trajectory (traj_cluster).
#'
add_clust_info <- function(f_graph, trajectories) {

  #Add info about clusters
  traj_cluster <- list()
  for(t in 1:length(trajectories)) {
    if (length(trajectories[[t]]$res) > 0) {
      traj_res_1 <- trajectories[[t]]$res[[1]]
      temp <- data.frame(node = as.vector(traj_res_1)) %>%
        left_join(f_graph$node_color, by = 'node')

      traj_cluster[[t]] <- data.frame(
        trajElmnts = paste(traj_res_1, collapse = " "),
        # clusterTraj = paste(unique(temp$cluster), collapse = ">")
        clusterTraj = paste(rle(as.vector(temp$cluster))$values, collapse = ">")
      )
    }
  }
  traj_clusters <- bind_rows(traj_cluster)
  traj_clusters
}
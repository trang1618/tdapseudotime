#' Simulated data from github.com/aridag/TDA_PSEUDOTIME.
#'
#' https://github.com/aridag/TDA_PSEUDOTIME/blob/master/MyDataSim.csv
#'
#' @format A data frame with 200 observations and 21 variables:
#' \code{id}, \code{covid_id}, \code{day}, \code{AlanineAminoTransferase},
#' \code{AspartateAminotransferase}, \code{Albumin},
#' \code{CardiacTroponinHighSensitivity}, \code{CRP}, \code{Creatinine},
#' \code{Ddimer}, \code{Ferritin}, \code{Fibrinogen}, \code{LDH}, \code{Lymphocyte},
#' \code{Neutrophil}, \code{Procalcitonin}, \code{PT}, \code{Bilirubin},
#' \code{WBC}, \code{first_date}, and \code{time}.
#'
#'
"sim_dat"

#' Matrix of the scaled lab values of \code{sim_dat}.
"scaled_lab_mat"

#' Data frame of centroids, each row represents the mean lab values of
#' all observations in each node of `f_graph`.
#' Output of the training workflow to establish the topology from fake data.
"centroids"

#' igraph object (from `graph.adjacency`)
#' output of the training workflow to establish the topology from fake data.
"f_graph"

#' List of trajectories computed from shortest paths.
#' output of the training workflow to establish the topology from fake data.
"out_trajectories"

#' Simulated data from github.com/aridag/TDA_PSEUDOTIME.
#'
#' https://github.com/aridag/TDA_PSEUDOTIME/blob/master/MyDataSim.csv
#'
#' @format A data frame with 23849 observations and 6 variables:
#' \code{siteid}, \code{patient_num}, \code{days_since_admission},
#' \code{concept_type}, \code{concept_code}, and \code{value}.
#'
#'
"sim_dat"

#' Processed data from sim_dat
#'
#' @format A data frame with 7207 observations and 19 variables:
#' \code{id}, \code{covid_id}, \code{time}, \code{AlanineAminoTransferase},
#' \code{AspartateAminotransferase}, \code{Albumin},
#' \code{CardiacTroponinHighSensitivity}, \code{CRP}, \code{Creatinine},
#' \code{Ddimer}, \code{Ferritin}, \code{Fibrinogen}, \code{LDH}, \code{Lymphocyte},
#' \code{Neutrophil}, \code{Procalcitonin}, \code{PT}, \code{Bilirubin},
#' \code{WBC}.
#'
#'
"processed_data"

#' Matrix of the scaled lab values of \code{sim_dat}.
"scaled_lab_mat"

#' Data frame of centroids, each row represents the mean lab values of
#' all observations in each node of `f_graph`.
#' Output of the training workflow to establish the topology from fake data.
"centroids"

#' Dataframe mapping node to cluster and corresponding color,
#' output of the training workflow to establish the topology from fake data.
"node_color"

#' List of trajectories computed from shortest paths.
#' output of the training workflow to establish the topology from fake data.
"out_trajectories"

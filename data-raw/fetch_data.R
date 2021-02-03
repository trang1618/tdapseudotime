set.seed(1618)
library(dplyr)

sim_dat <-
  read.csv(
    'https://raw.githubusercontent.com/aridag/TDA_PSEUDOTIME/7847a071c750ec1d80e179b8a5d6099b9d7b5573/PatientObservationsSim.csv'
  )
n_samps <- nrow(sim_dat)
patients <- unique(sim_dat$patient_num)
train_idx <- sample(patients, length(patients)*0.5, FALSE)
sim_dat %>%
  filter(patient_num %in% train_idx) %>%
  readr::write_csv('vignettes/train-data.csv')

sim_dat %>%
  filter(!patient_num %in% train_idx) %>%
  readr::write_csv('vignettes/test-data.csv')

###

processed_data <- widen_i2b2(sim_dat)
non_lab_value_names <- c('id', 'covid_id', 'CardiacTroponinHighSensitivity')
lab_value_names <- setdiff(names(processed_data), non_lab_value_names)
processed_data <- mutate_at(processed_data, dplyr::all_of(lab_value_names), as.numeric)

scaled_lab_mat <- processed_data[, lab_value_names] %>%
  mice(m = 1, maxit = 25, meth = 'pmm', seed = 500) %>% # imputation
  complete(1) %>%
  scale() %>% # prepare for cosine similarity calculation
  {.}

# scaled_lab_mat <- lapply(1:mm, function(x) complete(lab_values_mat, x)) %>%
#   do.call(rbind, .) %>%
#   as.matrix() %>%
#   scale()

usethis::use_data(sim_dat, scaled_lab_mat, processed_data, overwrite = TRUE)

# my_tda <- map_tda(scaled_lab_mat)
# my_graph <- make_tda_graph(my_tda, processed_data, 'time')

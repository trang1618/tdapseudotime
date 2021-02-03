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

usethis::use_data(sim_dat, overwrite = TRUE)

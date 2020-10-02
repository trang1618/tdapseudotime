set.seed(1618)
library(dplyr)

sim_dat <- read.csv('https://raw.githubusercontent.com/aridag/TDA_PSEUDOTIME/66dbdeb6239e0ff38947c3460a34374ea8895eda/MyDataSim.csv')
sim_dat <- sim_dat[sample(nrow(sim_dat), 400), ]
rownames(sim_dat) <- NULL

non_lab_value_names <- c('id', 'covid_id', 'day')
lab_value_names <- setdiff(names(sim_dat), non_lab_value_names)

sim_dat <- sim_dat %>%
  mutate(day = as.Date(day, format = '%Y-%m-%d')) %>%
  group_by(covid_id) %>%
  mutate(first_date = min(day)) %>%
  ungroup() %>%
  mutate(
    id = as.integer(id),
    time = as.numeric(day - first_date, units = 'days')) %>%
  arrange(covid_id, time) %>%
  mutate_at(dplyr::all_of(lab_value_names), as.numeric)

scaled_lab_mat <- scale(sim_dat[, lab_value_names])

usethis::use_data(sim_dat, scaled_lab_mat, overwrite = TRUE)
# cat(colnames(sim_dat), sep = '}, \\code{')

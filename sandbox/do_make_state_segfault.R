library(McMasterPandemic)
library(dplyr)

# doesn't fail
make_base_model(
  params = read_params("PHAC.csv"),
  state = make_state(params = read_params("PHAC.csv")),
  start_date = "2022-01-01",
  end_date = "2022-02-01",
  do_make_state = FALSE
)

# fails -- engages c++ make_state function and fails at the eigenvector step
make_base_model(
  params = read_params("PHAC.csv"),
  state = make_state(params = read_params("PHAC.csv")),
  start_date = "2022-01-01",
  end_date = "2022-02-01",
  do_make_state = TRUE
)

library(McMasterPandemic)
library(dplyr)

params <- read_params("ICU1.csv")
state <- make_state(params = params)

sir_state = (
  state[c("S", "Ia", "R")]
  %>% setNames(c("S", "I", "R"))
  %>% structure(epi_cat = c("S", "I", "R"), class = "state_pansim")
)
sir_params = (
  params[c('beta0', 'gamma_a')]
  %>% setNames(c("beta", "gamma"))
  %>% structure(class = "params_pansim")
)

M = matrix(0, nrow = 3, ncol = 3, dimnames = list(from = names(sir_state), to = names(sir_state)))


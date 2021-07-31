library("McMasterPandemic")

params1 <- read_params("ICU1.csv")
state1 <- make_state(params=params1)
M <- make_ratemat(params=params1, state=state1, sparse=TRUE)
s1A <- do_step(state1,params1, M, stoch_proc=FALSE)

s1A

library(TMB)
compile("TMB_do_step.cpp")
dyn.load(dynlib("TMB_do_step"))
library("McMasterPandemic")

params1 <- read_params("ICU1.csv")
state1 <- make_state(params=params1)
M <- make_ratemat(params=params1, state=state1, sparse=TRUE)
#s1A <- do_step(state1,params1, M, stoch_proc=FALSE)

dd <- MakeADFun(data = list(state = c(state1),
                            ratemat = M,
                            dt = 1,
                            do_hazard = TRUE,
                            stoch_proc = FALSE,
                            do_exponential = FALSE,
                            testwt_scale = "N"
                            ),
                parameters = list(params=c(params1)),
                DLL = "TMB_do_step")

new_state = dd$report()$state

new_state

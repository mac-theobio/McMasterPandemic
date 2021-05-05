library(McMasterPandemic)

params <- read_params("ICU1.csv")
paramsS <- update(params,c(proc_disp=0.1,obs_disp=100))
state <- make_state(params=params)
sdate <- "10-Feb-2020" ## arbitrary!
set.seed(101)
res1 <- run_sim(params,state,start_date=sdate,end_date="1-Jun-2020")
write.csv(res1, 'res1.csv', row.names=FALSE)

params <- read_params("ICU1.csv")
state <- make_state(params=params)
sdate <- "14-Feb-2020" ## arbitrary!
set.seed(101)
res1 <- run_sim(params,state,start_date=sdate,end_date="1-Jun-2020")
write.csv(res1, 'res2.csv', row.names=FALSE)

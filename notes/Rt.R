library(McMasterPandemic)
params <- fix_pars(read_params("ICU1.csv"))
params[["zeta"]] <- 1
rr1 <- run_sim(params=params, end_date="2020-Sep-1")
## incidence = foi*S
## beta(t) = incidence/(S*I) = foi/I
## why doesn't beta(t) -> 0?
## R(t) = beta(t)*S(0)
rr2 <- run_sim(params=update(params, zeta=8), end_date="2020-Sep-1")
with(rr,plot(date,foi*params[["N"]]/I,type="l",ylim=c(0,1)))
with(rr2,lines(date,foi*params[["N"]]/I, col=2))

with(rr,plot(date,incidence/S,type="l"))
with(rr,plot(date,incidence/S))

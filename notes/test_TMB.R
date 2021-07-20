library(TMB)
compile("TMB_do_step.cpp")
dyn.load(dynlib("TMB_do_step"))
library(McMasterPandemic)

p <- read_params("ICU1.csv")
s <- make_state(params = p)
s0 <- c(s)
init <- s0[1]
M <- make_ratemat(state = s, params = p, sparse=TRUE)
## note that attributes mess up MakeADFun - need to strip them with c()
## before passing to MakeADFun
print("***** Before calling MakeADFun ...")
dd <- MakeADFun(data = list(state = s0,
                            ratemat = M,
                            inf_ind = grep("I[a-z]", names(s)),
                            transm_ind = which(names(p) == "beta0"),
                            ## fragile! assumes same order as state
                            transm_wt_ind = grep("C[a-z]", names(p)),
                            foi_ind = c(which(rownames(M) == "S"),
                                        which(colnames(M) == "E"))

                            ),
                parameters = list(params=c(p)))

print("state in R = ")
print(s0)
print("ratemat in R = ")
print(M)
print("inf_ind in R= ")
print(grep("I[a-z]", names(s)))
print("transm_ind in R= ")
print(which(names(p) == "beta0"))
print("transm_wt_ind in R = ")
print(grep("C[a-z]", names(p)))
print("foi_ind in R= ")
print(c(which(rownames(M) == "S"), which(colnames(M) == "E")))

print("parameters in R = ")
print(list(params=c(p)))

print("***** Before calling cpp ...")
dd$fn(p)
print("***** After calling cpp ...")

identical(s0[1], init)
s == dd$report()$state ## new state

## it's probably possible to update data etc.
## by messing around with objects in the environment, e.g.
environment(dd$fn)$data$state
## ... rather than re-running MakeADFun() from scratch
## (we are departing a bit from the usual TMB use case,
## which is to optimize an objective function by calling
## the $fn and $gr elements many times with *constant*
## data (when we get to stage 2/3, do_step will be relegated
## to a function within a larger iteration/optimization process)

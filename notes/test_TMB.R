library(TMB)
compile("TMB_do_step.cpp")
dyn.load(dynlib("TMB_do_step"))
library(McMasterPandemic)

p <- read_params("ICU1.csv")
s <- make_state(params = p)
M <- make_ratemat(state = s, params = p, sparse=TRUE)

# Make the C++ function interface identical to do_step() in sim_funs.R
# NOTE: attributes mess up MakeADFun - need to strip them with c()
# before passing to MakeADFun

print("***** Before calling MakeADFun ...")
print(as.numeric(Sys.time())*1000000, digits=19)
dd <- MakeADFun(data = list(state = c(s),
                            ratemat = M,
                            dt = 1,
                            do_hazard = TRUE,
                            stoch_proc = FALSE,
                            do_exponential = FALSE,
                            testwt_scale = "Noooo!"
                            #inf_ind = grep("I[a-z]", names(s)),
                            #transm_ind = which(names(p) == "beta0"),
                            ## fragile! assumes same order as state
                            #transm_wt_ind = grep("C[a-z]", names(p)),
                            #foi_ind = c(which(rownames(M) == "S"),
                            #            which(colnames(M) == "E"))
                            ),
                parameters = list(params=c(p)))


print("***** Before calling cpp ...")
print(as.numeric(Sys.time())*1000000, digits=19)
#dd$fn(p)	# This doesn't call the C++ function
print("***** After calling cpp ...")
print(as.numeric(Sys.time())*1000000, digits=19)

s == dd$report()$state ## new state (this triggers the calling of the c++ function)

## it's probably possible to update data etc.
## by messing around with objects in the environment, e.g.
environment(dd$fn)$data$state
## ... rather than re-running MakeADFun() from scratch
## (we are departing a bit from the usual TMB use case,
## which is to optimize an objective function by calling
## the $fn and $gr elements many times with *constant*
## data (when we get to stage 2/3, do_step will be relegated
## to a function within a larger iteration/optimization process)

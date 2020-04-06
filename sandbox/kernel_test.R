devtools::load_all(".")
pars <- read_params(system.file("params","ICU1.csv",
                                package="McMasterPandemic"))

## FIXME: combination of numerical instability and epidemic phases
##   (pre-exponential, exponential, saturating, declining) makes this
##   a little sensitive to the length of the run
## The problem is likely related to allowing r to be negative
## Since it's a pandemic repo, we could just set the lower limit to zero for a while, I think.
## I am curious about BB's description of a "two-cycle", though
print(gg <- get_GI_moments(pars))
nt <- gg[["Gbar"]]*10
kk <- transKernel(pars, do_hazard=FALSE, steps=nt)$foi
print(kernelMoments(kk))

pars[["E0"]] <- 0.001
sim <- run_sim_range(params=pars, nt=gg[["Gbar"]]*5)
rI <- diff(log(sim$Is))
print(max(tail(rI,10)))


devtools::load_all(".")
pars <- read_params(system.file("params","ICU1.csv",
                                package="McMasterPandemic"))

print(get_GI_moments(pars))
try(print(kernelMoments(transKernel(pars, do_hazard=FALSE, steps=500)$foi)))
    
sim <- run_sim(pars, end_date="01-Jul-2020")
rI <- diff(log(sim$Is))
print(max(rI[-1:-10]))


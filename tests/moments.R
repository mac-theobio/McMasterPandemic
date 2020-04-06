library(McMasterPandemic)

## here the two different Gbar computations are close
## (about 0.3% difference)
params <- read_params(system.file("params","ICU1.csv",
                                  package="McMasterPandemic"))
get_GI_moments(params)
get_kernel_moments(params)


## this example comes from a calibration to r=0.1914 ..., Gbar=6
## calibration was done with get_GI_moments, kernel moments don't match?
## which should we trust?
## numerical issues?  make dt smaller?
ccS <- readRDS(system.file("testdata","ccS.rds",
                           package="McMasterPandemic"))
get_GI_moments(ccS$param)
get_kernel_moments(ccS$param)

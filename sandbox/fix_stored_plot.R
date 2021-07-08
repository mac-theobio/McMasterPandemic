library(McMasterPandemic)
source("../tests/testthat/utils.R")

mod <- readRDS("ontario_vaccine.RDS")
fix_mod <- fix_stored(mod$fit)

plot(fix_mod,data=mod$fitdat)

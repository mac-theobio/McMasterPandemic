library(McMasterPandemic)

mod <- readRDS("ontario_vaccine.RDS")

plot(mod$fit,data=mod$fitdat)

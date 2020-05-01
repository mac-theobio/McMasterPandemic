library(mvbutils)
library(McMasterPandemic)
ff <- foodweb(where=2,rprune="forecast|calibrate|run_",plotting=FALSE, ancestors=FALSE)
plot(ff,lwd=2)

library(mvbutils)
library(McMasterPandemic)
ff <- foodweb(where=2,rprune="forecast|calibrate|run_")
plot(ff,lwd=2)

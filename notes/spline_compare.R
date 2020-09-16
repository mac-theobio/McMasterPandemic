library(McMasterPandemic)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(plot(ff, data=dat))

coef(ff, "fitted")

print(mod_ns)



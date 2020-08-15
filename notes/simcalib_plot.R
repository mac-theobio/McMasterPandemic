library(ggplot2);theme_set(theme_bw())
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

gg <- (ggplot(dat, aes(x=date, color=var))
	+ geom_line(aes(y=value))
	+ geom_point(aes(y=data))
	+ facet_grid(testing_intensity~keepvars, scale="free")
	+ scale_y_log10()
)

print(gg %+% filter(dat,opt_testify==TRUE) + ggtitle("Opt testify"))
print(gg %+% filter(dat, opt_testify == FALSE) + ggtitle("No opt testify"))





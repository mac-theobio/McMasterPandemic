library(ggplot2);theme_set(theme_bw())
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(head(dat))

gg <- (ggplot(dat, aes(x=date, color=var))
	+ geom_line(aes(y=value))
	+ geom_point(aes(y=data),alpha=0.1)
	+ facet_wrap(~testing_type, scale="free")
	+ scale_y_log10()
)


print(gg)




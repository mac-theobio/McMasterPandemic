library(ggplot2);theme_set(theme_bw())
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(head(dat))

print(tail(dat))


gg <- (ggplot(dat, aes(x=date, color=var))
	+ geom_line(aes(y=value))
	+ geom_point(aes(y=data),alpha=0.4)
	+ facet_wrap(~var, scale="free", ncol=2)
#	+ scale_y_log10(limits=c(1,NA))
	+ scale_y_log10()
	+ ylab("Daily counts")
	+ ggtitle("Logistic testing")
)

print(gg)



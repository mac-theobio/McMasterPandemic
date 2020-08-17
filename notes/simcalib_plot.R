library(ggplot2);theme_set(theme_bw())
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(head(dat))

gg <- (ggplot(dat, aes(x=date, color=var))
	+ geom_line(aes(y=value))
	+ geom_point(aes(y=data),alpha=0.4)
	+ facet_wrap(~testing_type, scale="free", ncol=1)
	+ scale_y_log10(limits=c(1,NA))
	+ ylab("Daily counts")
)

print(gg)

ggpostest <- (ggplot(dat %>% filter(var=="postest"), aes(x=date,color=testing_type))
	+ geom_line(aes(y=value))
	+ geom_point(aes(y=data),alpha=0.4)
	+ scale_y_log10(limits=c(1,NA))
	+ scale_color_manual(values=c("black","red","blue"))
	+ ylab("Daily counts")
	+ ggtitle("Positive test")
)

print(ggpostest)


library(ggplot2);theme_set(theme_bw())
library(dplyr)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(simdat)

gg <- (ggplot(simdat, aes(x=date,y=med, color=type))
	+ geom_line()
	+ geom_point(aes(shape=type))
   + geom_ribbon(aes(ymin=lwr,ymax=upr,fill=type), alpha=0.3)
	+ scale_color_manual(values=c("black","red"))
	+ scale_fill_manual(values=c("black","red"))
	+ facet_wrap(~var,scale="free",ncol=1)
	+ scale_y_log10()
	+ ylab("Daily count")
)

print(gg)

print(gg %+% (simdat %>% filter(date < as.Date("2020-04-15"))))



library(ggplot2);theme_set(theme_bw())
library(dplyr)
library(cowplot)

source("makestuff/makeRfuns.R")
commandFiles()

gg <- (ggplot(simcombo, aes(x=date,y=med, color=type, shape=var, group=type))
	+ geom_point()
	+ scale_color_manual(values=c("black","red","blue"))
	+ scale_fill_manual(values=c("black","red","blue"))
	+ scale_y_log10()
	+ ylab("Daily count")
	+ facet_wrap(~params, ncol=3)
	+ theme(legend.position = "bottom")
)


print(gg)

print(gg %+% (simcombo %>% filter(var != "incidence")))


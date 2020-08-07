library(ggplot2);theme_set(theme_bw())
library(dplyr)
library(cowplot)

source("makestuff/makeRfuns.R")
commandFiles()

gg <- (ggplot(simdat, aes(x=date,y=value, color=var))
	+ geom_point()
	+ scale_color_manual(values=c("black","red","blue","orange"))
	+ scale_y_log10()
	+ ylab("Daily count")
	+ facet_grid(W_asymp~iso_p)
	+ theme(legend.position = "bottom")
)


print(gg)


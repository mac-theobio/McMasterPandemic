library(ggplot2);theme_set(theme_bw())
library(dplyr)
library(cowplot)

source("makestuff/makeRfuns.R")
commandFiles()

gg <- (ggplot(simdat)
	+ aes(x=date,y=value)
	+ geom_line()
	+ scale_color_manual(values=c("black","red","blue","orange"))
	+ scale_y_log10(limits=c(1, NA))
	+ ylab("Daily count")
	+ theme(legend.position = "bottom")
)

print(gg + facet_grid(W_asymp~iso_p) + aes(color=var))
print(gg + facet_grid(W_asymp~var) + aes(color=iso_p))
print(gg + facet_grid(iso_p~var) + aes(color=W_asymp))

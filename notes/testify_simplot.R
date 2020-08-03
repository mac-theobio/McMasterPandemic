library(ggplot2);theme_set(theme_bw())
library(dplyr)
library(cowplot)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(simdat)

gg <- (ggplot(simdat, aes(x=date,y=med, color=type))
	+ geom_line()
	+ geom_point(aes(shape=type))
   + geom_ribbon(aes(ymin=lwr,ymax=upr,fill=type), alpha=0.3)
	+ scale_color_manual(values=c("black","red","blue"))
	+ scale_fill_manual(values=c("black","red","blue"))
	+ facet_wrap(~var,scale="free",nrow=1)
	+ scale_y_log10()
	+ ylab("Daily count")
	+ ggtitle("Default MacPan testify")
)

ggdefault <- gg + theme(legend.position="none")

ggwts <- (gg 
	%+% simdatwts 
	+ theme(legend.position="bottom") 
	+ ggtitle("Different testing weights")
)

ggPos <- (gg
	%+% simdatPos
	+ theme(legend.position="none")
	+ ggtitle("Different positivity")
)

ggcombo <- plot_grid(ggdefault,ggPos,ggwts,ncol=1)

print(ggcombo)


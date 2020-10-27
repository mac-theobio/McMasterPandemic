library(ggplot2);theme_set(theme_bw())
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()
print(ensembles_dat)

fulldat <- ensembles_dat %>% filter(trim == "full") 
shortdat <- ensembles_dat %>% filter(trim == "short")


gg <- (ggplot(data=filter(shortdat,var == "report")
	, aes(date))
	+ geom_point(aes(y=value),size=0.1)
	+ geom_line(aes(y=med),size=0.1)
	+ geom_ribbon(aes(ymin=lwr, ymax=upr),alpha=0.3) 
	+ facet_grid(spline_df~spline_pen)
)

print(gg + ggtitle("short fit reported cases") + scale_y_log10())
print(gg %+% filter(shortdat,var == "Rt") + ggtitle("short fit Rt"))

#print(gg %+% filter(fulldat,var == "report")
#	+ ggtitle("full fit reported cases")
#	+ scale_y_log10()
#	)

#print(gg %+% filter(fulldat,var == "Rt")
#	+ ggtitle("full fit Rt")
#	+ scale_y_continuous(limit=c(0,4), oob=scales::squish)
#)






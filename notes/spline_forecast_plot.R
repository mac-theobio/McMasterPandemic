library(ggplot2);theme_set(theme_bw())
library(dplyr)

source("makestuff/makeRfuns.R") ## Will eventually be a package
commandEnvironments() ## Read in any environments specified as dependencies
makeGraphics()

print(head(forecastdat))


gg <- (ggplot(data=forecastdat, aes(x=date))
#	+ scale_y_log10()
	+ geom_point(aes(y=value,color=type,alpha=type),size=0.1)
	+ scale_color_manual(values=c("red","black"))
	+ scale_alpha_manual(values=c(0.3,1))
	+ geom_ribbon(aes(ymin=lwr,ymax=upr,fill=type),alpha=0.3)
	+ facet_wrap(~seed,scale="free")
	+ scale_fill_manual(values=c("red","black"))
)

print(gg)

print(gg + scale_y_log10())

print(gg %+% filter(forecastdat, date>=as.Date("2020-05-30")))

## To save as .rda
## saveVars(y) ## OR
## saveEnviroment()

## To save as .rds
## rdsSave(x)

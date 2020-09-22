library(ggplot2);theme_set(theme_bw())
library(dplyr)

source("makestuff/makeRfuns.R") ## Will eventually be a package
commandEnvironments() ## Read in any environments specified as dependencies
makeGraphics()

print(head(forecastdat))

forecastdat <- forecastdat %>% mutate(seed=ifelse(is.na(seed),101,seed))
simdat <- forecastdat %>% filter(type == "sim")
truedat <- forecastdat %>% filter(type == "true")

gg <- (ggplot(data=simdat, aes(x=date))
#	+ geom_line(aes(y=value),color="red")
	+ geom_point(data=truedat, aes(x=date,y=value),color="black",size=0.1)
	+ scale_y_log10()
	+ geom_ribbon(aes(ymin=lwr,ymax=upr,fill=type),alpha=0.5)
	+ facet_wrap(~seed,scale="free")
)

print(gg)


## To save as .rda
## saveVars(y) ## OR
## saveEnviroment()

## To save as .rds
## rdsSave(x)

library(ggplot2);theme_set(theme_bw()+theme(panel.spacing=grid::unit(0,"lines")))
library(dplyr)
library(tidyr)
library(cowplot)

callArgs <- "combine_testing.Rout combine_testing.R testwt_N.random.clean_testing.rda testwt_N.symptomatic.clean_testing.rda testwt_N.focus.clean_testing.rda"


source("makestuff/makeRfuns.R")
commandEnvironments()
if (!interactive()) makeGraphics()

simdat1 <- (bind_rows(
	getEnvironment("random")[["simdat1"]]
	, getEnvironment("symptomatic")[["simdat1"]]
	, getEnvironment("focus")[["simdat1"]]
	)
	%>% mutate(testing_intensity = factor(testing_intensity,levels=c("2 Test per 1000 per day","10 Test per 1000 per day"))
		, testing_type = factor(testing_type, levels = c("Random", "Symptomatic", "Proactive")
		)
	)
)

simdat2 <- (bind_rows(
	getEnvironment("random")[["simdat2"]]
	, getEnvironment("symptomatic")[["simdat2"]]
	, getEnvironment("focus")[["simdat2"]]
	)
	%>% mutate(testing_intensity = factor(testing_intensity
		, levels=c("2 Test per 1000 per day", "10 Test per 1000 per day"))
		, testing_type = factor(testing_type, levels = c("Random", "Symptomatic", "Proactive")
		)
	)
)

ggall <- (ggplot(simdat1
	, aes(x=date,y=value,color=testing_type,linetype=isolation)
	)
 	+ scale_colour_manual(values=c("black","red","blue"))
	+ geom_line()
	+ facet_grid(strength~testing_intensity, scale="free")
	+ theme(legend.position = "bottom")
)

print(ggall 
	%+% (simdat1 
		%>% filter(varname == "Cumulative Incidence")
		)
	+ ggtitle("Cumulative Incidence")
)


ggall2 <- (ggplot(simdat2, aes(x=date,y=value,color=testing_type, linetype = isolation))
 	+ scale_colour_manual(values=c("black","red","blue"))
	+ geom_line()
	+ facet_grid(strength2~testing_intensity2, scale="free")
	+ theme(legend.position = "bottom")
)

print(ggall2 
	%+% (simdat2 %>% filter(varname == "Cumulative Incidence")
	)
	+ ggtitle("Cumulative Incidence")
)


print(ggall
	%+% (simdat1
		%>% filter(varname == "Positivity")
	)
	+ ggtitle("Positivity")
)

print(ggall2
	%+% (simdat2
		%>% filter(varname == "Positivity")
	)
	+ ggtitle("Positivity")
)

saveVars(simdat1, simdat2)

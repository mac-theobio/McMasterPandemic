library(ggplot2);theme_set(theme_bw()+theme(panel.spacing=grid::unit(0,"lines")))
library(dplyr)
library(tidyr)
library(cowplot)
library(shellpipes)

callArgs <- "testwt_N.random.simreportplot.Rout testify_simreportplot.R testwt_N.random_sims.rda"
source("makestuff/makeRfuns.R")
commandEnvironments()
if (!interactive()) makeGraphics()


if(grepl("random",targetname()))testing_wts <- "Random"
if(grepl("symptomatic",targetname()))testing_wts <- "Symptomatic"
if(grepl("focus",targetname()))testing_wts <- "Proactive"


varnames <- data.frame(var = c("total_test", "positivity", "postest", "CumIncidence", "incidence")
	, varname = c("Daily Test", "Positivity", "Positive tests\nper day", "Cumulative Incidence", "Incidence")
)

simdat <- (simdat 
	%>% left_join(.,varnames)
	%>% rename(isolation = iso_t)
	%>% rename(strength = R0)
	%>% mutate(isolation = as.character(isolation)
		, testing_intensity = paste0(testing_intensity*1000, " Test per 1000 per day")
		, strength = ifelse(strength == 2, "fast", "slow")
	)
	%>% ungroup()
)

simdat1 <- (simdat 
	%>% filter(testing_type == "constant")
	%>% mutate(testing_type = testing_wts)
)

print(dim(simdat1))

ggall <- (ggplot(simdat1
	, aes(x=date,y=value,linetype=isolation)
	)
# + scale_colour_manual(values=c("red","blue","black"))
	+ geom_line()
	+ facet_grid(strength~testing_intensity, scale="free")
)

print(ggall 
	%+% (simdat1 
		%>% filter(var == "CumIncidence")
		)
	+ ggtitle(testing_wts)
)

simdat2 <- (simdat
	%>% filter(testing_intensity == "2 Test per 1000 per day")
	%>% filter(strength == "fast")
	%>% mutate(testing_intensity2 = ifelse(testing_type %in% c("mix","increaseT"), "Increasing Testing","Constant Testing")
		, strength2 = ifelse(testing_type %in% c("mix","reduceR"), "Decreasing R", "Constant R")
	)
	%>% select(-testing_type)
	%>% mutate(testing_type = testing_wts)
#	%>% distinct()
)

print(dim(simdat2))

ggall2 <- (ggplot(simdat2, aes(x=date,y=value, linetype = isolation))
	+ geom_line()
	+ facet_grid(strength2~testing_intensity2, scale="free")
	+ ggtitle(testing_wts)
)

print(ggall2)

saveVars(simdat1, simdat2)

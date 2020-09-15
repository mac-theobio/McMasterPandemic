library(McMasterPandemic)
library(tidyverse)


source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

params <- read_params("ICU1.csv")


dat <- (sims 
	%>% select(date, report, death)
	%>% gather(key = "var", value="value", -date)
)

ff <- calibrate_comb(params
	





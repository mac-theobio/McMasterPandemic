library(McMasterPandemic)
library(tidyverse)


source("makestuff/makeRfuns.R") ## Will eventually be a package
commandEnvironments() ## Read in any environments specified as dependencies
makeGraphics()

flist <- list.files(path="cachestuff/",pattern="_calib[.]")

print(flist)

plot_fits <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	print(plot(modlist$fit, data=modlist$fitdat))
	print(modlist$fitdat)
}

sapply(flist,plot_fits)



library(McMasterPandemic)

source("makestuff/makeRfuns.R") ## Will eventually be a package
commandEnvironments() ## Read in any environments specified as dependencies
## makeGraphics()


flist <- list.files(path="cachestuff/",pattern="recalib")

print(flist)



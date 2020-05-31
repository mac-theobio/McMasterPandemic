library(McMasterPandemic)
library(splines)
library(dplyr)
library(parallel)

L <- load(".ont_keep.RData")
print(unique(ont_all_sub$var))
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")
params <- fix_pars(read_params("ICU1.csv"))
dat <- ont_noICU
datH <- dat %>% filter(var == "H")
datreport <- dat %>% filter(var == "report")



Honly <- do.call(calibrate_comb,c(nlist(params, data=datH, use_spline=TRUE, use_phenomhet=TRUE,DE_optim=FALSE)))
reportonly <- do.call(calibrate_comb,c(nlist(params, data=datH, use_spline=TRUE, use_phenomhet=TRUE,DE_optim=FALSE)))

plot(Honly)
plot(reportonly)

#rdsave(reportonly,Honly)


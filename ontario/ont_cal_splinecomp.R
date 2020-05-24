library(McMasterPandemic)
library(splines)
library(dplyr)
library(parallel)

L <- load(".ont_keep.RData")
print(unique(ont_all_sub$var))
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")
params <- fix_pars(read_params("ICU1.csv"))

inputs <- expand.grid(spline_df=3:7,
                      spline_sd_pen=c(NA,0.001,0.01,0.1),
                      knot_quantile_var=c(NA,"H","report"))

res_list <- mclapply(seq(nrow(inputs)),
                     function(i) {
                         cat(i,"\n")
                         do.call(calibrate_comb,
                                 c(nlist(params
                                       , debug_plot=FALSE
                                       , data=ont_noICU)
                                 , inputs[i,]))
                     },
                     mc.cores=5)

#rdsave(res_list)


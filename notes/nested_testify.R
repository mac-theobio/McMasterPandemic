library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(parallel)
library(furrr)
library(future.batchtools)

testing_scale <- 0.2
max_intensity <- 0.002
max_intensity <- 0.5


callArgs <- "basic.nested_testify.Rout nested_testify.R batchtools.rda testify_funs.rda basic.rda sims.csv"

source("makestuff/makeRfuns.R")
print(commandEnvironments())
makeGraphics()

fn <- if (interactive()) "PHAC_testify.csv" else matchFile(".csv$")
pars <- (read_params(fn)
    %>% fix_pars(target=c(R0=R0, Gbar=Gbar))
    %>% update(N=pop)
)

## different combinations
testing_intensity <- c(0.001)
keep_vars <- c("postest/H/death")
opt_testify <- c(FALSE)
constant_testing <- c(TRUE,FALSE)

comboframe <- expand.grid(testing_intensity=testing_intensity
	, keep_vars = keep_vars
	, opt_testify = opt_testify
	, constant_testing = constant_testing
)

datevec <- as.Date("2020-01-01"):as.Date("2020-10-01")
testdat <- data.frame(Date = as.Date(datevec)
	# , intensity = plogis(seq(-1,1,length.out = length(datevec)),scale=testing_scale)*max_intensity
	, intensity = seq(0.0001,0.001,length.out = length(datevec))
)
plot(testdat)

dd <- (simtestify(p=update_pars(comboframe[1,]),testdat)
   %>% transmute(date,H,death,postest)
   %>% gather(key="var",value="value",-date)
)

ggplot(dd,aes(x=date,y=value,color=var))+geom_point()+scale_y_log10(limits=c(1,NA))


sim_and_calibrate <- function(y,testdat){
	x <- comboframe[y,]
   if(x[4]){
      testdat$intensity <- 1
   }
	pp <- update_pars(x)
	simdat <- simtestify(pp,testdat)
	calib_mod <- calibrate_sim(dd=simdat, pars=pp, p=x, testdat)
#	calib_mod <- NULL
	res_list <- list(fit=calib_mod,params=pp, data=simdat)
	saveRDS(object=res_list, file=paste0("./cachestuff/simcalib.",y,".RDS"))
	return(res_list)
}

batch_setup()

res_list <- future_map(seq(nrow(comboframe)),function(x)sim_and_calibrate(x,testdat))

## interactive playing around stuff

if (FALSE) {
    debug(simulate_testify_sim)
    debug(calibrate_sim)
    sim_and_calibrate(1)
    do.call(calibrate_comb, c(nlist(params = pars, use_DEoptim = FALSE, 
                                    use_spline = FALSE,
                                    debug=TRUE,
                                    debug_plot=TRUE,
                                    data = dat2, opt_pars = opt_pars, sim_args = sim_args, 
                                    maxit = 1000)))
    ## params.log_beta0                params.log_E0 
    ## -0.6579844                    1.7804940 
    ## params.log_testing_intensity          log_nb_disp.postest 
    ## -9.1067633                    2.0989440

    ## beta0                E0 testing_intensity 
    ## -0.5062236         1.6094379        -6.9077553 


}

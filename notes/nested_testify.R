library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(parallel)
library(furrr)
library(future.batchtools)

callArgs <- "testwt_N.nested_testify.Rout nested_testify.R batchtools.rda testify_funs.rda testwt_N.rda sims.csv"

source("makestuff/makeRfuns.R")
print(commandEnvironments())
# makeGraphics()

## prevent breakage of old input files when we add new parameters
default_vals <- list(stoch_obs=FALSE,
                     obs_disp=1,
                     testing_intensity=c(0.001),
							testing_type = c("constant","linear","logistic"),
                     keep_vars=c("postest/H/death"),
                     opt_testify=FALSE
                     )

## if it wasn't already specified, set it from the default value

for (nm in names(default_vals)) {
    if (!exists(nm)) assign(nm,default_vals[[nm]])
}
    

pars <- (read_params(matchFile(".csv$"))
    %>% fix_pars(target=c(R0=R0, Gbar=Gbar))
    %>% update(N=pop)
)

## If we are not optimizing for it, we need to set min/starting to true
pars[["testing_intensity"]] <- testing_intensity 

## single dispersion parameter for *all* observed vars ...
if (stoch_obs) {
    pars <- update(pars, obs_disp=obs_disp)
}

## different combinations

comboframe <- expand.grid(keep_vars = keep_vars
	, opt_testify = opt_testify
	, testing_type = testing_type
)


print(comboframe)

datevec <- as.Date(start):as.Date(end)
testdat <- data.frame(Date = as.Date(datevec)
	, intensity = testing_intensity
)

plot(testdat)

dd<- (simtestify(p=pars, testdat)
   %>% select(date,H,death,postest)
   %>% gather(key="var",value="value",-date)
   %>% mutate(type="time varying")
)
print(ggplot(dd,aes(x=date,y=value,color=type))
   + geom_line()
   + facet_wrap(~var,ncol=2)
   + scale_y_log10(limits=c(1,NA))
)

sim_and_calibrate <- function(y,testdat,debug_plot){
	x <- comboframe[y,]
	if(x$testing_type == "linear"){
		testdat$intensity <- seq(min_testing,max_testing,length.out =nrow(testdat))
	}
	if(x$testing_type == "logistic"){
      testdat$intensity <- plogis(seq(qlogis(min_testing/max_testing),qlogis(0.99),length.out = nrow(testdat)))*max_testing
	}
	simdat <- simtestify(pars,testdat)
	calib_mod <- calibrate_sim(dd=simdat, pars=pars, p=x, testdat,
                                   debug_plot=debug_plot, debug=TRUE, debug_hist=TRUE)
#	calib_mod <- NULL
	res_list <- list(fit=calib_mod,params=pars, data=simdat)
        ## BMB: is this OK or is copying/moving stuff into cachestuff supposed to be done make-ily?
	saveRDS(object=res_list, file=paste0("./cachestuff/simcalib.",y,".RDS"))
	return(res_list)
}
 
# res <- sim_and_calibrate(1,testdat, debug_plot=TRUE)
# 
# ## plot parameter histories
# hh <- (attr(res$fit,"debug_hist")
#     %>% as_tibble()
#     %>% mutate(n=seq(nrow(.)))
#     %>% pivot_longer(-n)
#     %>% mutate_at("name", ~forcats::fct_inorder(factor(.)))
# )
# ggplot(hh,aes(n,value,colour=name))+facet_wrap(~name,scale="free_y") + geom_line()
## stop here
# quit()
batch_setup()

print(comboframe)
print(seq(nrow(comboframe)))

res_list <- future_map(seq(nrow(comboframe)),function(x) sim_and_calibrate(x,testdat,debug_plot=FALSE))

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


saveVars(comboframe)

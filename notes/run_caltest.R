library(McMasterPandemic)
library(tidyverse)
library(anytime)
library(bbmle)
library(parallel)

sim_cali <- function(seed, do_plots=FALSE,use_DEoptim=TRUE) {
	cat(seed,"\n")
   set.seed(seed)
   simdat <- forecast_sim(p=true_pars
   	, opt_pars = opt_pars
		, base_params = params
		, start_date=start_date
   	, end_date=end_date
   	, stoch = c(obs = TRUE, proc=FALSE)
   # 	, time_args = list(mob_value=mob_vals   ## mob_value can't be same length as TS
   # 		, mob_startdate=mob_start
   # 		)
   #    , sim_fun=run_sim_mobility
	)
    ## plot(simdat, log=TRUE)
   simdat <- (simdat
   %>% filter(var %in% c("report","hosp"))
	# %>% filter(!is.na(value)) 
   %>% mutate(value = ifelse(is.na(value),0,value))
   )
   
   simdat2 <- simdat %>% filter(between(date, cut_start, cut_end))
   
    ## print(params)
    ## print(opt_pars)
    g1 <- do.call(calibrate_comb,
                c(nlist(params=params
                        , debug_plot=FALSE
                        ## whatever data you need to have modified
                        , data=simdat2
                        , opt_pars=opt_pars
                )))

    print(bbmle::coef(g1$mle2))

    return(list(data=simdat2,fit=g1,fulldata=simdat))
}

## sim_cali(1)
res <- mclapply(seq(nsim), sim_cali,mc.cores = 1)
# res <- lapply(seq(nsim), sim_cali)
# res2 <- lapply(seq(nsim),function(x)sim_cali(seed=x,use_DEoptim=TRUE))

print(res)
print(res2)
# rdsave(res,res2,true_pars)

library(McMasterPandemic)
library(tidyverse)
library(anytime)
library(bbmle)
library(parallel)

if (!exists("use_DEoptim")) use_DEoptim <- FALSE

sim_cali <- function(seed, do_plots=FALSE,use_DEoptim=FALSE) {
	cat(seed,"\n")
   set.seed(seed)
   simdat <- forecast_sim(p=true_pars
   	, opt_pars = opt_pars
		, base_params = params
		, start_date=start_date
   	, end_date=end_date
   	, stoch = c(obs = TRUE, proc=FALSE)
   	, time_args = list(mob_value=mob_vals   ## mob_value can't be same length as TS
   		, mob_startdate=mob_start
   		)
      , sim_fun=run_sim_mobility
	)
    ## plot(simdat, log=TRUE)
   simdat <- (simdat
   %>% filter(var %in% c("report"))
	%>% filter(!is.na(value)) 
   )
   
   simdat <- simdat %>% filter(between(date, cut_start, cut_end))
   
   if (do_plots) plot(ggplot(simdat,aes(date,value))+geom_point()+scale_y_log10())

    ## print(params)
    ## print(opt_pars)
    g1 <- calibrate(data=simdat, base_params=params
    	, start_date = cut_start
      , opt_pars = opt_pars
      , debug_plot=FALSE
      ## , debug=FALSE
      , use_DEoptim = use_DEoptim
      , DE_cores = 1 ## don't fight with mclapply
      ## , mle2_args=list(browse_obj=TRUE)
      , time_args = list(mob_value=mob_vals
      	, mob_startdate=mob_start)
      , sim_fun=run_sim_mobility
      )

    print(bbmle::coef(g1$mle2))

    pp <- predict(g1)

    res_dat <- data.frame(bbmle::confint(g1$mle2, method="quad", level=0.95)
                        , estimate = bbmle::coef(g1$mle2)
                        , seed = seed
                        , pars = names(g1$mle2@coef)
                          )
    print(true_pars)
    print(res_dat)
    return(list(simdat=simdat,fit=g1,pars=res_dat,pred=pp, fullsim=simdat))
}

## sim_cali(1)
#res <- mclapply(seq(nsim), sim_cali)
res <- mclapply(seq(nsim), sim_cali)
res2 <- lapply(seq(nsim),function(x)sim_cali(seed=x,use_DEoptim=TRUE))

print(res)
print(res2)
# rdsave(res,res2,true_pars)

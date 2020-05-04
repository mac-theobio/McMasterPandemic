library(McMasterPandemic)
library(tidyverse)
library(anytime)
library(bbmle)
library(parallel)

if (!exists("use_DEoptim")) use_DEoptim <- FALSE

sim_cali <- function(seed, do_plots=FALSE) {
	cat(seed,"\n")
   set.seed(seed)
   simdat <- forecast_sim(p=true_pars
   	, opt_pars = opt_pars
		, base_params = params
		, start_date=start_date
   	, end_date=end_date
   	, break_dates = bd
   	, stoch = c(obs = TRUE, proc=FALSE)
	)
    ## plot(simdat, log=TRUE)
   simdat <- (simdat
   %>% filter(var %in% c("report"))
	%>% filter(!is.na(value)) 
	%>% mutate(value = round(value))
   )
   
   simdat <- simdat %>% filter(between(date, cut_start, cut_end))
   
   if (do_plots) plot(ggplot(simdat,aes(date,value))+geom_point()+scale_y_log10())

    ## print(params)
    ## print(opt_pars)
    g1 <- calibrate(data=simdat, base_params=params
                  , start_date = cut_start
                  , opt_pars = opt_pars
                  , break_dates = bd
                  ## , debug_plot=TRUE
                  ## , debug=FALSE
                  , use_DEoptim = use_DEoptim
                  , DE_cores = 1 ## don't fight with mclapply
                    ## , mle2_args=list(browse_obj=TRUE)
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
res <- mclapply(seq(nsim), sim_cali)

print(res)

# rdsave(res,true_pars)

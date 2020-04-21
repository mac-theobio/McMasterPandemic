library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)

## setup 

params <- fix_pars(read_params("ICU1.csv"))
params[["obs_disp"]] <- 5

start_date <- anydate("2020-01-01")
end_date <- anydate("2020-04-01")

cutoff_start <- anydate("2020-02-10")
cutoff_end <- anydate("2020-03-10")

bd <- "2020-02-22"
bd <- NULL

opt_pars <- list(
   params=c(log_E0=4, log_beta0=-1)
   , log_rel_beta0 = rep(-1, length(bd))
   , log_nb_disp=0
)

params[["N"]] <- 1e7

nsim <- 20

sim_cali <- function(x){
   
   set.seed(x)
   sim1break <- run_sim(params
      , start_date=start_date
      , end_date=end_date
      , stoch = c(obs = TRUE, proc=FALSE)
   )
   
   simdat <- pivot(condense(sim1break))
   simdat <- simdat %>% filter(var %in% c("report"))
   
   ## Need to round value because of negative binomial fit
   dd <- (simdat 
          %>% filter(!is.na(value)) 
          %>% mutate(value = round(value))
          %>% filter(between(date,cutoff_start, cutoff_end))
   )
   
   g1 <- calibrate(data=dd, base_params=params
                   , opt_pars = opt_pars
                   , break_dates = bd
   )
   
   res_dat <- data.frame(bbmle::confint(g1$mle2, method="quad", level=0.95)
                         , estimate = bbmle::coef(g1$mle2)
                         , seed = x
                         , pars = names(g1$mle2@coef)
   )
   
   return(res_dat)
}

ll <- lapply(1:nsim,function(x)sim_cali(x))

truedf <- data.frame(pars = c("params.log_E0","params.log_beta0","log_rel_beta0","log_nb_disp")
                     , trueval = c(log(params[["E0"]]), log(params[["beta0"]]), log(rel_break1), log(params[["obs_disp"]]))
)

simdf <- (bind_rows(ll)
          %>% left_join(.,truedf)  
          %>% rowwise()
          %>% transmute(pars
                        , trueval
                        , estimate
                        , lwr = X2.5..
                        , upr = X97.5..
                        , inCI = between(trueval, lwr, upr)
                        , seed
          )
          %>% ungroup()
          %>% group_by(pars)
          %>% mutate(ind = rank(estimate))
)

ggmilli <- (ggplot(simdf, aes(x=ind,y=estimate,color=inCI))
            + geom_point()
            + geom_pointrange(aes(ymin=lwr,ymax=upr))
            + geom_hline(aes(yintercept = trueval))
            + facet_wrap(~pars,scale="free")
            + theme_bw()
            + scale_color_manual(values=c("red","blue","black"))
)
print(ggmilli)

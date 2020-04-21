library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)

L <- load("run_caltest.RData")
## setup 

params <- fix_pars(read_params("ICU1.csv"))
params[["obs_disp"]] <- 5

start_date <- anydate("2020-01-01")
end_date <- anydate("2020-04-01")

# cutoff_start <- anydate("2020-02-10")
# cutoff_end <- anydate("2020-03-10")

cutoff_start <- anydate("2020-01-01")
cutoff_end <- anydate("2020-04-01")

break1 <- "2020-02-22"
bd <- anydate(c(break1))
rel_break1 <- 0.2

opt_pars <- list(
   ## these params go to run_sim
   params=c(log_E0=4, log_beta0=-1)
   , log_rel_beta0 = rep(-1, length(bd))
   , log_nb_disp=0
)
params[["N"]] <- 1e7

truedf <- data.frame(pars = c("params.log_E0","params.log_beta0","log_rel_beta0","log_nb_disp")
                   , trueval = log(c(params[c("E0","beta0")]
                                    , rel_break1
                                    , params[["obs_disp"]])))


names(res) <- seq_along(res) ## seeds
simvals <- (map_dfr(res,~left_join(rename(.$simdat,sim=value),
                                   rename(.$pred,pred=value),
                                   by=c("date","var")),.id="seed")
    %>% select(-vtype)
    %>% pivot_longer(names_to="type",cols=c("sim","pred"))
)

ggplot(simvals, aes(date,value,lty=type)) +
    geom_line() +
    facet_wrap(~seed) +
    scale_y_log10() +
    geom_vline(xintercept=bd,colour="red")

simdf <- (map_dfr(res,pluck,"pars")
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

simwide <- (simdf
    %>% select(pars,estimate,seed)
    %>% pivot_wider(names_from="pars",values_from="estimate")
)

pairs(simwide[,-1],gap=FALSE)

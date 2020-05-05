library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)
library(bbmle)
library(DEoptim)

## DEoptim calibrated object is not working with the pipeline

# load("run_DEoptim_breaks.RData")
# res <- res[c(1,3,5,7,9,11,13,15,17,19)]
## Old working example

# load("run_caltest_nobreak.RData")

truedf <- data.frame(pars = names(true_pars)
	, trueval = true_pars
)

rescombo <- c(res,res2)

names(rescombo) <- seq_along(rescombo) ## seeds


likvals <- suppressWarnings(
    map_dfr(rescombo,~tibble(NLL=-bbmle::logLik(.$fit$mle2)),.id="seed")
)

print(likvals)

# seed_order <- likvals %>% arrange(NLL) %>% pull(seed)

simvals <- (map_dfr(rescombo,~left_join(rename(.$simdat,sim=value),
                                   rename(.$pred,pred=value),
                                   by=c("date","var")),.id="seed")
#    %>% mutate_at("seed", factor,levels=seed_order)
    %>% select(-vtype)
    %>% pivot_longer(names_to="type",cols=c("sim","pred"))
)

gg1 <- ggplot(simvals, aes(date,value,lty=type)) +
    geom_line() +
    facet_wrap(~seed)

## print(gg1)
# print(gg1 + scale_y_log10())

simdf <- (map_dfr(rescombo,pluck,"pars")
    %>% left_join(.,truedf)
    %>% rename(lwr="X2.5..",upr="X97.5..")
    %>% mutate(inCI =
                   case_when(is.na(lwr) ~ "?",
                             lwr < trueval & trueval < upr ~ "yes",
                             TRUE ~ "no"
                             )
               )
)

print(simdf)
simdf$type <- rep(c("NM","DEoptim"),each=4*length(res))

#simrank <- (simdf
#    %>% filter(pars=="params.log_beta0")
#    %>% mutate(ind =rank(estimate))
#    %>% select(seed, ind)
#)
#simdf <- full_join(simdf, simrank, by="seed")

ggmilli <- (ggplot(simdf, aes(x=seed,y=estimate,color=type))
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

## all 'converged' except fit 3
map_dbl(res, ~.$fit$mle2@details$convergence)
## local conv?

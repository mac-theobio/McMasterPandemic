library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)
library(bbmle)
library(DEoptim)

trueNLL <- function(x){
	m <- McMasterPandemic:::mle_fun
	f <- x[["fit"]]$forecast_args
	f <- f[!names(f) == "fixed_pars"]
	dd <- (x[["simdat"]]
   	%>% filter(!is.na(value)) 
    	%>% mutate(value = round(value))
	)
	NLL <- do.call(m,c(list(true_pars,data=dd),f))
	return(NLL)
}


truedf <- data.frame(pars = names(true_pars)
	, trueval = true_pars
)

rescombo <- c(res,res2)

names(rescombo) <- seq_along(rescombo) ## seeds

NLL <- map_dbl(rescombo,~-logLik(.$fit$mle2))
NLLtrue <- map_dbl(rescombo,~trueNLL(.))

LLdat <- data.frame(seed = c(1:length(res),1:length(res))
	, type = rep(c("Nelder-Mead","DEoptim"),each=length(res))
	, NLL = NLL
	, NLLtrue = NLLtrue
	, NLL_diff = NLL - NLLtrue
)
	
	
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

ggLL <- (ggplot(LLdat, aes(x=seed, y=NLL_diff,color=type))
	+ geom_point()
	+ geom_hline(yintercept = 0)
	+ scale_color_manual(values = c("red","blue"))
	+ theme_bw()
)

print(ggLL)

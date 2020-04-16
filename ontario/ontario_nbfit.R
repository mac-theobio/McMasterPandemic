library(dplyr)
library(glmmTMB)
library(broom.mixed)
library(ggplot2); theme_set(theme_bw())

fitfun <- function(d,key) {
    fit <- glmmTMB(value~poly(vday,2,raw=TRUE),
                    family=nbinom2,
                   data=d)
    ## assume that non-pos-def result is caused by nbinom ...
    if (is.na(AIC(fit))) {
        fit <- update(fit, family=poisson)
    }
    return(fit)
}
ont_recent_ntv <- filter(ont_recent_nt, var != "Ventilator")
ont_nb_fit <- setNames(dplyr::group_map(group_by(ont_recent_ntv,var), fitfun),
                  unique(ont_recent_ntv$var))
fix_term <- function(x) {
  gsub("poly(vday, 2, raw = TRUE)","Poly" , x, fixed=TRUE)
}
ont_nb_tab <- (purrr::map_dfr(ont_nb_fit,tidy,conf.int=TRUE,.id="var")
    %>% mutate_at("term",fix_term)
    %>% mutate(type=case_when(term=="(Intercept)" ~ "int", 
                              grepl("1$", term) ~ "linear",
                              TRUE ~ "quad")
             , var=factor(var,levels=rev(
                                  c("newConfirmations","Hospitalization","ICU","newDeaths")
                              ))
               )
    %>% select(-c(effect,component,statistic,p.value))
)

gg3 <- (ggplot(ont_nb_tab,aes(y=var,x=estimate,xmin=conf.low,xmax=conf.high))
    + geom_pointrange()
    + facet_wrap(~type,scale="free",ncol=1)
    + geom_vline(xintercept=0,lty=2)
    + labs(y="")
)

print(gg3)

# rdsave("ont_nb_fit","ont_nb_tab")

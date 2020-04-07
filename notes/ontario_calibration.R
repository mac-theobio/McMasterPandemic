library(tidyverse)
library(glmmTMB)
library(nlme)
library(broom.mixed)
url <- "https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/clean.Rout.csv"
dd <- read_csv(url)

## process with knitr::spin() ...

##'  Filter/date by time since first reported case, per province ...

##' Warm up with a quadratic fit (changing scale in ggplot link function in GLMs) to reported incidence in all provinces
##+ first_plot
print(gg0 <- ggplot(dd,aes(Date,incidence,colour=Province))
    + geom_line(size=0.5,alpha=0.3)
    + geom_point()+scale_y_log10()
    + geom_smooth(method="lm",
                  formula=y~poly(x,2),
                  se=FALSE)
)

## easiest way to plot mixed model outputs with CIs?
## sjPlot, broom.mixed, ... ?

## data exploration and possible comparison of IDEA/IHME/etc. approaches
## (i.e. purely phenomenological)
## fit linear, quad mixed? models
## logit, Richards?
## epigrowthfit or mle2 or ... ?
## plot predictions and CIs

ont_dd <- (dd
    %>% filter(Province=="ON")
    %>% select(Date,Hospitalization,ICU,ventilator,incidence,deceased)
    %>% pivot_longer(-Date,names_to="var")
)

print(gg1 <- ggplot(ont_dd,aes(Date,value,colour=var))
    + geom_point()
    + scale_y_log10()
    + geom_smooth(method="lm",
                  formula=y~poly(x,2))
)
##' * Natural for ICU and ventilator to be parallel (no real lags here)
##' * recent nonlinearity/flattening in cases is too recent/sharp to affect the estimated slope much
##' * very weird that the ICU/vent curves are flattening before the other two

##' Picture is cleaner (and flattening is more apparent) if we focus on more recent data:
##'

## find first useful day
min_day <- function(day,value) {
    good <- !is.na(value) & value>0
    if (all(good)) min(day) else min(day[good])
}

ont_recent <- (filter(ont_dd,Date>=as.Date("2020-03-15"))
    %>% mutate(day=as.numeric(Date-min(Date)))
    %>% group_by(var)
    %>% mutate(vday=day-min_day(day,value))
    %>% ungroup()
)
    
print(gg0 %+% ont_recent)

print(gg2 <- ggplot(ont_recent,aes(vday,value,colour=var))
    + geom_point()
    + scale_y_log10()
    + geom_smooth(method="lm",
                  formula=y~poly(x,2))
    + scale_x_continuous(limits=c(-1,NA))
)

fit1 <- glmmTMB(value~var-1 + var:vday,
        family=nbinom2,
        dispformula=~var-1,
        data=ont_recent)

## FIXME: check vs nbinom1 ?

##' use raw polynomials to parameterize in terms of initial slope
##' rather than average slope
fit1Q <- update(fit1, . ~ -1 + var + var:poly(vday,2,raw=TRUE))
## false convergence warning?

fix_term <- function(x) {
    gsub(":poly(vday, 2, raw = TRUE)","Poly",
         gsub("var","",x),
         fixed=TRUE)
}

t0 <- (tidy(fit1Q,conf.int=TRUE)
    %>% mutate_at("term",fix_term)
    %>% mutate(type=ifelse(grepl("[[:alpha:]]$",term),"int",
                       ifelse(grepl("1$",term),"linear","quad")))
    %>% mutate(var=gsub("Poly[0-9]","",term),
               ## reverse levels to get top-to-bottom
               var=factor(var,levels=rev(c("incidence","Hospitalization",
                                       "ICU","ventilator","deceased"))))
    %>% select(-c(effect,component,statistic,p.value))
)
(gg2 <- ggplot(t0,aes(y=var,x=estimate,xmin=conf.low,xmax=conf.high))
    + geom_pointrange()
    + facet_wrap(~type,scale="free",ncol=1)
    + geom_vline(xintercept=0,lty=2)
)

##' * intercepts for variables other than incidence describe sensitivity (how small a non-zero value is actually reported?
##' * incidence and hosp slopes and curvature agree fairly well
##' * no 'clear' evidence of slowing incidence/hospitalization (¿ flattening in incidence could be obscured by increasing testing/ too recent or sharp to be well-modeled by a quadratic ?)
##' * 
##' * don't know why death is slower  and ICU/ventilator are faster
##' * deaths conf intervals might be overly narrow because these are cumulative values

##' calibrate to hospitalization

vv <- vcov(fit1Q)$cond
hosp_ind <- grep("Hospitalization",rownames(vv))
vv <- vv
## implied curves of R_t?

## calibrate to reports


## dates of control???

## run sims


### ISSUES


## 1. non-constant r (use epigrowthfit and feed in observed relative r over time (i.e. derivative of fitted curve) as timevars?
## 2. non-constant testing (ignore this for now, maybe we can do a brute-force correction where we subtract r(testing) from r(reports) - this assumes that tests are not expanding to incorporate a different subset of people
## 3. use H, ICU, D to get r then use reports to get time-varying beta0?

## suppose we have GLM(M)NB fit to something
## * sample from the distribution of slopes and intercepts
## for each sample of (i0, r):
##    * adjust the parameters to get the right beta0/r
##    * maybe pick a random value of G?
##    * brute-force calibrate initial conditions based on those parameters and the i0 we picked
## 
##    * run the simulation to forecast

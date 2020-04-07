library(tidyverse)
library(glmmTMB)
library(nlme)
url <- "https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/clean.Rout.csv"
dd <- read_csv(url)

## Filter/date by time since first reported case, per province ...

## %>% mutate(t_abs=as.numeric(Date-min(Date)))
## %>% group_by(Province)
## %>% mutate(t_loc=as.numeric(Date-min(Date)

## quadratic fit (changing scale in ggplot confuses
##  link function in GLMs)
(ggplot(dd,aes(Date,incidence,colour=Province))
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

## implied curves of R_t?

## calibrate to reports

## dates of control???

## run sims



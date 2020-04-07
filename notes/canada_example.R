library(tidyverse)
library(glmmTMB)
library(nlme)
url <- "https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/clean.Rout.csv"
dd <- read_csv(url)

## Filter/date by time since first reported case, per province ...

## %>% mutate(t_abs=as.numeric(Date-min(Date)))
## %>% group_by(Province)
## %>% mutate(t_loc=as.numeric(Date-min(Date)

ggplot(dd,aes(Date,incidence,colour=Province))+
    geom_point()+scale_y_log10()+
    geom_smooth(method="glm",
                formula=y~poly(x,2),
                se=FALSE)

## fit linear, quad mixed? models
## logit, Richards?


## calibrate to reports

## dates of control???

## run sims



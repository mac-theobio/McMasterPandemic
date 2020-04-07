library(tidyverse)
library(glmmTMB)
url <- "https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/clean.Rout.csv"


ggplot(dd,aes(Date,incidence,colour=Province))+
    geom_point()+scale_y_log10()+
    geom_smooth(method="glm",
                formula=y~poly(x,2),
                se=FALSE)
                
## fit linear, quad mixed? models
## logit, Richards?

dd <- read_csv(url)

## calibrate to reports

## dates of control???

## run sims



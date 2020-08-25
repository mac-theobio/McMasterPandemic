## Using epigrowthfit to estimate growth rate with incidence/hospitalization/ICU/ventilator data from Ontario, Canada
library(epigrowthfit)
library(tidyverse)
library(glmmTMB)
library(nlme)
library(broom.mixed)

url <- "https://raw.githubusercontent.com/wzmli/COVID19-Canada/master/clean.Rout.csv"
dd <- read_csv(url)

start_date <- "2020-03-08"

## is epigrowthfit going to complain about NAs? .... ans: YES!
## Too many NA for hospitalization, let's try ICU and ventilators. Do I even want to fit it to prevalance?
## need to fill in NA; fill in ICU and ventilator for 2020-03-29
## ICU = 81, ventilator = 54

ICUfill <- 81
Ventfill <- 54

ont_dd <- (dd
  %>% filter(Province=="ON")
  %>% filter(Date >= as.Date(start_date))
  %>% select(Date,Hospitalization,ICU,ventilator,new_report=incidence,deceased)
  %>% mutate(time = 1:nrow(.)
    , ICU = ifelse(Date == as.Date("2020-03-29"), ICUfill,ICU)
    , ventilator = ifelse(Date == as.Date("2020-03-29"), Ventfill,ventilator)
    , new_hosp = diff(c(NA,Hospitalization))  ## Do I even want the new_XX?
    , new_icu = diff(c(NA,ICU))
    , new_vent = diff(c(NA,ventilator))
    )
  %>% filter(Date >= as.Date(start_date))
)


reportfit <- epigrowthfit(data=ont_dd
  , deaths = ont_dd$new_report
  , optimizer = "nlminb"
  , verbose = TRUE
  , model = "logistic"
  , drop_mle2_call = FALSE
  , optCtrl = list(eval.max=1e9, iter.max=1e9)
)

reportfit@mle2@details$convergence

plot(reportfit)


## Logistic converged and gave a slightly larger r compared to exp (0.187 vs. 0.1725)




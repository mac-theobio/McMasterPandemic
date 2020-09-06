library(McMasterPandemic)
library(tidyverse)

pp <- read_params("ICU1.csv") %>% fix_pars()
r <- run_sim(params=pp, end_date=as.Date("2020-06-15"))

## X is the hospital accumulator (cumulative total of hospitalizations)
plot(r,keep_states=c("X","H","hosp"), log=TRUE)

r2 <- r %>% select(date,X,H,hosp)
head(r2)
identical(r$hosp,c(NA,diff(r$X)))

## now with testify
ppt <- read_params("PHAC_testify.csv") %>% fix_pars()
rt <- run_sim(params=ppt, end_date=as.Date("2020-06-15"),
             ratemat_args=list(testify=TRUE))

## X is the hospital accumulator (cumulative total of hospitalizations)
plot(rt,keep_states=c("X","H","hosp"), log=TRUE)

r2t <- rt %>% select(date,X,H,hosp)
identical(rt$hosp,c(NA,diff(rt$X)))


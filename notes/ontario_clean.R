library(readr)
library(dplyr)
library(tidyr)

## FIXME: should cache results
url <- "https://wzmli.github.io/COVID19-Canada/git_push/clean.Rout.csv"
dd <- read_csv(url)
PHOurl <- "http://wzmli.github.io/COVID19-Canada/PHO.csv"
ddPHO <- read_csv(PHOurl)
ont_dd <- (dd
  %>% filter(Province=="ON")
  %>% select(Date,Hospitalization,ICU,Ventilator,deceased,newConfirmations)
    %>% pivot_longer(-Date,names_to="var")
    %>% setNames(tolower(names(.)))
)
min_day <- function(day,value) {
    good <- !is.na(value) & value>0
    if (all(good)) min(day) else min(day[good])
}

start_date <- "2020-03-15"

ont_recent <- (ont_dd
    %>% filter(date>=as.Date(start_date))
    %>% mutate(day=as.numeric(date-min(date)))
    %>% group_by(var)
    %>% mutate(vday=day-min_day(day,value))
    %>% ungroup()
)

save("ont_recent","start_date","ont_dd", file="ontario_clean.RData")

## use rdsave ... ?

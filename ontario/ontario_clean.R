library(readr)
library(dplyr)
library(tidyr)

## FIXME: should cache results?
url <- "https://wzmli.github.io/COVID19-Canada/git_push/clean.Rout.csv"
dd <- read_csv(url)
## filter(dd,Date==as.Date("2020-04-14"),Province=="ON")
PHOurl <- "http://wzmli.github.io/COVID19-Canada/PHO.csv"
ddPHO <- read_csv(PHOurl)
ont_all <- (dd
    %>% filter(Province=="ON")
    %>% select(Date,Hospitalization,ICU,Ventilator,deceased,newConfirmations,newTests)
    %>% mutate(newDeaths=c(NA,diff(deceased)),
               ## ON hosp includes ICU, our model compartment is just acute care
               Hospitalization=Hospitalization-ICU)
    %>% select(-deceased)
    %>% pivot_longer(-Date,names_to="var")
    %>% setNames(tolower(names(.)))
)
min_day <- function(day,value) {
    good <- !is.na(value) & value>0
    if (all(good)) min(day) else min(day[good])
}

start_date <- "2020-03-15"
## start_date <- "2020-03-01"

ont_recent <- (ont_all
    %>% filter(date>=as.Date(start_date))
    %>% mutate(day=as.numeric(date-min(date)))
    %>% group_by(var)
    %>% mutate(vday=day-min_day(day,value))
    %>% ungroup()
)

ont_recent_nt <- filter(ont_recent,var != "newTests")
ont_all_nt <- filter(ont_all,var != "newTests")

save("ont_all", file=sprintf("ONdata_%s.rda",
                              format(Sys.time(),"%Y%b%d")))
# rdsave("ont_all", "ont_recent", "ont_all_nt", "ont_recent_nt", "start_date")


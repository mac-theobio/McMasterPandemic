library(McMasterPandemic)
library(readr)
library(dplyr)
library(tidyr)
library(zoo)

source("makestuff/makeRfuns.R")
commandEnvironments()

tsdat_url <- "https://wzmli.github.io/COVID19-Canada/git_push/clean.Rout.csv"
	
tsdat <- read_csv(tsdat_url)

### Clean ts data
Ontario_dat <- (tsdat
	%>% filter(Province=="ON")
   %>% select(Province,Date,Hospitalization,ICU,Ventilator,deceased,newConfirmations,newTests)
	%>% mutate(newDeaths=c(NA,diff(deceased))
   	## ON hosp includes ICU, our model compartment is just acute care
   	, Hospitalization=Hospitalization-ICU)
   %>% select(-deceased)
   %>% pivot_longer(names_to="var",-c(Date,Province))
   %>% setNames(tolower(names(.)))
	%>% ungroup()
)

## translate variable names to internally used values
## drop unused variables
keep_vars <- c("H","ICU","death","report","newTests")

## Maybe keep reports only for simplicity

keep_vars <- c("report")

clean_tsdata <- (Ontario_dat
    %>% mutate_at("var",trans_state_vars)
    %>% filter(var %in% keep_vars)
)

date_vec <- as.Date(min(clean_tsdata$date):max(clean_tsdata$date))

date_df <- data.frame(date = rep(date_vec,1)
   , var = rep(c("report"),each=length(date_vec))
   )

ont_dat <- (left_join(date_df,clean_tsdata))

print(ont_dat)

saveVars(ont_dat)


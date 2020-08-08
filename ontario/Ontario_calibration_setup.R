library(McMasterPandemic)
library(tidyverse)
library(zoo)

source("makestuff/makeRfuns.R")
## Don't worry about this when coding interactively/sourcing
commandEnvironments()


## Section 1: Read Data Sources

### Get time series data from MLi's data repo; mobility data from Apple and Google.

tsdat_url <- "https://wzmli.github.io/COVID19-Canada/git_push/clean.Rout.csv"
google_url <- "https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv"
apple_url <- "https://raw.githubusercontent.com/ActiveConclusion/COVID19_mobility/master/apple_reports/applemobilitytrends.csv"
	
tsdat <- read_csv(tsdat_url)
apple <- read_csv(apple_url)
google <- read_csv(google_url)

## Section 2: Clean data
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

clean_tsdata <- (Ontario_dat
    %>% mutate_at("var",trans_state_vars)
    %>% filter(var %in% keep_vars)
)

### Clean Mobility data

clean_google <- (google
    %>% filter(country_region == "Canada", sub_region_1 == "Ontario")
    %>% select(date, contains("baseline"))
    %>% pivot_longer(names_to="type", values_to="value", -c(date))
    %>% mutate_at("date", as.Date)
    %>% mutate_at("type", str_remove, "\\_percent.*")
    %>% mutate_at("value", ~./100+1)
)

clean_apple <- (apple
    %>% filter(region  == "Ontario")
    %>% rename(subregion="sub-region")
    %>% select(-c(geo_type,country,alternative_name,subregion,region))
    %>% pivot_longer(names_to="date", values_to="value",
                     -c(transportation_type))
    ## names_transform=list(date=as.Date))
    %>% mutate_at("date", as.Date)
    %>% mutate_at("value", ~./100)
    %>% rename(type="transportation_type")
)

clean_mobility <- (bind_rows(google=clean_google,apple=clean_apple,.id="source")
    %>% mutate(tvec=as.numeric(date-min(date,na.rm=TRUE)))
    %>% filter(type %in% c("retail_and_recreation","workplaces","driving"))
    %>% dplyr::select(date,value)
    %>% group_by(date)
    %>% summarise_at("value",mean,na.rm=TRUE)
    %>% na.omit()
    %>% rename(rel_activity = value)
#   %>% mutate_at("rel_activity", ~pmin(., 1))  ## cap at 100% (? should we ?)
	 %>% ungroup()
)

## Section 3 McMaster calibration setup

### filter ts data to calibrate

calibrate_data <- (clean_tsdata
   %>% filter(var %in% c("H","report","death"))
   %>% mutate(value = ifelse(is.na(value),0,value))
)

date_vec <- as.Date(min(calibrate_data$date):max(calibrate_data$date))

date_df <- data.frame(date = rep(date_vec,3)
   , var = rep(c("H","report","death"),each=length(date_vec))
)

## Why does it take Date and not date? Does it matter?
test_df <- data.frame(date = date_vec
	, var = rep("newTests",length(date_vec))
)

test_data_fill <- (left_join(test_df, clean_tsdata)
	%>% transmute(Date = date
		, intensity = ifelse(is.na(value),0,value)
	)
)

print(test_data_fill)


calibrate_data_fill <- (left_join(date_df,calibrate_data)
   %>% mutate(value = ifelse(is.na(value),0,value))
)

saveVars(calibrate_data_fill, test_data_fill, clean_mobility, ext="rda")


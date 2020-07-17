library(McMasterPandemic)
library(tidyverse)

## July 16th


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
Ontario_dat <- (dat
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
   %>% mutate_at("rel_activity", ~pmin(., 1))  ## cap at 100% (? should we ?)
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

calibrate_data_fill <- (left_join(date_df,calibrate_data)
   %>% mutate(value = ifelse(is.na(value),0,value))
)

### MacPan setup
params <- fix_pars(read_params("ICU1.csv")
   , target=c(Gbar=6)
   , pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
)

params[["N"]] <- 14.57e6 ## Population of Ontario (2019)

### parameters we are trying to estimate
opt_pars <- list(params=c(log_E0=2, log_beta0=-1, logit_c_prop=-1, logit_mu = -1, logit_phi1 = -1),log_nb_disp=c(report=1, death=1,H=1))


## Section 4: Calibrate 

Ontario_fit <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
      , data=calibrate_data_fill
      , mob_data = clean_mobility
      , opt_pars = opt_pars
      , use_DEoptim = TRUE
		, DE_cores = 2
		, use_phenomhet = TRUE
		, use_mobility = TRUE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
#		, use_spline , spline_df , spline_days ## Not using splines right now
	)
	)
)

save.image(Ontario_fit, calibrate_data_fill, clean_mobility, file = "basic.rda")


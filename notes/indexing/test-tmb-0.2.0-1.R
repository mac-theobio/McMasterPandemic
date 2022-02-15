# negative log likelihood

library(ggplot2)
library(McMasterPandemic)
library(lubridate)
library(tidyr)
library(dplyr)
set_spec_version("0.2.0", '../../inst/tmb')
tmb_mode()

load_model_from_file = TRUE
options(MP_force_symm_vcov = TRUE)

params1 <- read_params("PHAC.csv")
params1[c("N", "phi1")] <- c(42507, 0.98)
state1 <- make_state(params=params1)

# start and end dates
sdate <- as.Date("2021-11-01")
edate <- as.Date("2022-01-19")
initial_date = as.Date("2021-08-03")
start_date_offset = as.integer(sdate - initial_date)

# read and process data
covid_data <- ("../../sandbox/yukon/report_data_yukon_h_and_i.csv"
               %>% read.csv
               %>% mutate(date = as.Date(date))
               %>% filter(date >= ymd(20210803))
               %>% filter(between(as.Date(date), sdate, edate))

               # report -- new reported cases on that day
               # hosp -- new hospital admissions on that day
               %>% select(date, report, hosp)

               %>% pivot_longer(names_to = "var", -date)
               %>% mutate(value=round(value))
)

head(covid_data, n=12)

# establish schedule of time variation of parameters
params_timevar = data.frame(
  Date = ymd(
    # estimate a new transmission rate on
    # these dates (i'm no expert but these
    # seemed to "work")
    20211115, # nov 15 beta0
    20211215, # dec 15 beta0
    20211215, # dec 15 mu
    20220101 # jan 01 beta0
  ),
  Symbol = c("beta0", "beta0", 'mu', 'beta0'),
  Value = NA,
  Type = "abs"
)

mm = (make_base_model(
    params = params1,
    start_date = sdate - start_date_offset,
    end_date = edate,
    params_timevar = params_timevar,
    data = covid_data
  )
  %>% update_tmb_indices
)
sim_args = list(flexmodel = mm)
mm$tmb_indices$observed

obj_fun = tmb_fun(mm)
obj_fun$env$data  # items available to be read in to tmb using DATA_* macros

## reprex for calibrate() opt_pars vs time_args issue

## based off of calibration demo from getting_started vignette

library(McMasterPandemic) ## should be installed from master branch

## get synthetic data with a single change in beta0
params1 <- read_params("ICU1.csv")
state1 <- make_state(params=params1)

sdate <- "2020-02-10"
edate <- "2020-06-01"

set.seed(101)
params1obs <- update(params1, obs_disp=200)

params_timevar <- data.frame(
  Date = "2020-04-01",
  Symbol = "beta0",
  Relative_value = 0.5
)

res1obs <- run_sim(params1obs, state1,
                   start_date=sdate,
                   end_date=edate,
                   params_timevar = params_timevar,
                   stoch=c(obs=TRUE, proc=FALSE),
                   step_args=list(do_hazard = TRUE))

report_data <- (res1obs
                %>% mutate(value=round(report), var="report")
                %>% select(date, value, var)
                %>% na.omit()
)

## initial conditions for calibrating the single change in beta0
opt_pars <- list(time_params = c(0.1))

## calibration, setting value in *both* opt_pars and time_args$params_timevar
fitted.mod <- calibrate(
  data = report_data
  , start_date = sdate
  , base_params = params1obs
  , opt_pars = opt_pars
  ## "oops" accidentally trying to fix the value we're also trying to calibrate
  ## who will win??
  , time_args = list(params_timevar = params_timevar)
)

print(paste0("true relative beta0: ", params_timevar$Relative_value))

print(paste0("initial calibration value for relative beta0: ", opt_pars$time_params))

print(paste0("estimated beta0 from calibrate: ", coef(fitted.mod, "fitted")))

## calibration, forgetting to pass time_args$params_timevar
fitted.mod <- calibrate(
  data = report_data
  , start_date = sdate
  , base_params = params1obs
  , opt_pars = opt_pars
)

print(paste0("true relative beta0: ", params_timevar$Relative_value))

print(paste0("initial calibration value for relative beta0: ", opt_pars$time_params))

print(paste0("estimated beta0 from calibrate: ", coef(fitted.mod, "fitted")))

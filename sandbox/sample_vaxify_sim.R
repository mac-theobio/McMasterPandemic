devtools::load_all()
library(dplyr)
library(ggplot2)

## sim settings
start_date <- "2020-02-01"
end_date <- "2020-06-01"
time_pars <- data.frame(Date=c("2020-03-01"),
                        Symbol=c("vax_doses_per_day"),
                        Relative_value = c(1.5))
# time_pars <- data.frame(Date=c("2020-03-15", "2020-05-10"),
#                          Symbol=c("beta0"),
#                          Adjustment_type = c("relative", "relative"),
#                          Adjustment_value=c(0.7, 1.4))

## params
## set up vax params
vax_doses_per_day <- 1e5

## initialize params
base_params <- update(read_params("PHAC.csv")
                       , N = 1.4e7
)

# vax_params <- update(expand_params_vax(base_params,
#                                         vax_doses_per_day = vax_doses_per_day),
#                       obs_disp = 10)
vax_params <- update(expand_params_vax(base_params,
                                       vax_doses_per_day = vax_doses_per_day))

## initialize states
vax_state <- make_state(params = vax_params)

# res_vax <- run_sim(vax_params, vax_state,
#                    start_date = start_date, end_date = end_date,
#                    params_timevar = time_pars,
#                    stoch = c(obs = TRUE, proc = FALSE),
#                    step_args = list(do_hazard = TRUE),
#                    condense_args = list(keep_all = TRUE))
res_vax <- run_sim(vax_params, vax_state,
                   start_date = start_date, end_date = end_date,
                   params_timevar = time_pars,
                   stoch = c(obs = TRUE, proc = FALSE),
                   step_args = list(do_hazard = TRUE),
                   condense_args = list(keep_all = TRUE))

ggplot(get_doses_per_day(res_vax), aes(x = date, y = total_doses_per_day)) + geom_point()

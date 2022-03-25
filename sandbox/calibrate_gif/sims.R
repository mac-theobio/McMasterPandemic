library(McMasterPandemic)
library(dplyr)
library(animation)

reset_spec_version()
tmb_mode()

params <- read_params("ICU1.csv")
start_date = "2021-05-10"
end_date = "2021-12-10"
# tv_dat <- data.frame(
#   Date = c("2021-06-18", "2021-07-01", "2021-07-25"),
#   Symbol = c("beta0", "beta0", "beta0"),
#   Value = c(0.1, 0.5, 0.9),
#   Type = c("rel_orig", "rel_orig", "rel_orig")
# )

model <- make_base_model(params,
                         state = make_state(params = params),
                         start_date = start_date, end_date = end_date,
#                         params_timevar = tv_dat,
                         do_hazard = TRUE)

r_sim <- run_sim(
  params = c(model$params, obs_disp = 5),
  state = model$state,
  start_date = start_date,
  end_date = end_date, stoch = c(obs = TRUE, proc = FALSE),
#  params_timevar = tv_dat,
  step_args = list(do_hazard = TRUE),
  condense = TRUE
)

report_data <- (r_sim
                %>% mutate(value = round(report), var = "report")
                %>% select(date, value, var)
                %>% na.omit()
)
if(FALSE) plot(report_data$value, type = "l")

params[["beta0"]] = 0.3

ani.options(interval = 0.1, nmax = 100)

saveGIF(
  fitted_tmb <- calibrate(
    data = report_data,
    time_args = list(), #list(params_timevar = tv_dat),
    base_params = params,
    debug_plot = TRUE,
    mle2_control = list(maxit = 15),
    opt_pars = list(params = c(beta0 = params[["beta0"]])),
    sim_args = list(
      step_args = list(do_hazard = TRUE),
      use_flex = TRUE,
      flexmodel = model
    )
  )
)

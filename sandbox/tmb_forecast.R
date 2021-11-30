library(McMasterPandemic)
options(macpan_pfun_method = "grep")
library(dplyr)
## FIXME: need to import dplyr bind_rows()
## FIXME: can we do.call(rbind, ...)) instead???
## FIXME: consider poorman (won't work for pivoting though)
## we also depend on %>% (don't do it!) if we
## require R >= 4.0.0 we can just use |> instead

options(macpan_pfun_method = "grep")
  options(MP_use_state_rounding = FALSE)
  options(MP_vax_make_state_with_hazard = FALSE)

  start_date <- "2021-02-01"
  end_date <- "2021-09-01"

  params_timevar = data.frame(
    Date = as.Date(c("2021-04-20", "2021-06-20", "2021-08-20", "2021-05-03", "2021-07-08")),
    Symbol = c("beta0", "beta0", "beta0", "Ca", "Ca"),
    Value = c(0.5, 0.1, 2.2, 0.2, 0.01),
    Type = c("rel_orig")
  ) %>% arrange(Symbol, Date)

  base_params <- read_params("PHAC.csv")
  vax_params <- expand_params_vax(
    params = base_params,
    vax_doses_per_day = 1e4,
    model_type = "twodose"
  )
  model_params <- expand_params_variant(
    vax_params,
    variant_prop = 0.5,
    variant_advantage = 1.5,
    variant_vax_efficacy_dose1 = 0.3,
    variant_vax_efficacy_dose2 = 0.8
  ) %>% expand_params_S0(1 - 1e-5)

  state = make_state(params = model_params)

  model = make_vaccination_model(
    params = model_params,
    state = state,
    start_date = start_date, end_date = end_date,
    do_hazard = TRUE,
    do_approx_hazard = FALSE,
    do_make_state = TRUE,
    max_iters_eig_pow_meth = 100,
    tol_eig_pow_meth = 1e-6,
    params_timevar = params_timevar,
    do_variant = TRUE)

  time_wrap(
    tmb_sim <- run_sim(
      params = model_params,
      state = state,
      start_date = start_date, end_date = end_date,
      params_timevar = params_timevar,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      use_flex = TRUE,
      flexmodel = model
    )
   ,
    r_sim <- run_sim(
      params = model_params,
      state = state,
      start_date = start_date, end_date = end_date,
      params_timevar = params_timevar,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      use_flex = FALSE
    )
  )


## FIXME: why do we need to specify start-date and
## end-date in two different places??

## FIXME: can use_flex be '(!is.null(flexmodel))'
##  instead?

## FIXME: SCW; fix case for params_timevar = NULL

## FIXME: assign class to result properly

tmb_sim <- forecast_sim(
    p = NULL,
    opt_pars = list(),
    base_params = model_params,
    start_date = start_date,
    end_date = end_date,
    time_args = list(params_timevar = params_timevar),
    sim_args = list(
        step_args = list(do_hazard = TRUE),
        use_flex = TRUE,
        condense = FALSE,
        flexmodel = model)
)


test_sim <- function(use_flex = TRUE, obs_disp = NA) {
  arg_list <- list(
      p = NULL,
      opt_pars = list(),
      base_params = model_params,
      start_date = start_date,
      end_date = end_date,
      time_args = list(params_timevar = params_timevar),
      sim_args = list(
          step_args = list(do_hazard = TRUE),
          use_flex = use_flex,
          condense = FALSE,
          flexmodel = model)
  )
  if (!is.na(obs_disp)) {
    arg_list$base_params <- c(arg_list$base_params,
                              list(obs_disp = obs_disp))

    arg_list <- c(arg_list, list(stoch = c(obs=TRUE)))
  }
  do.call(forecast_sim, arg_list)
}

t1 <- test_sim() ## no stoch
t3 <- test_sim(use_flex=FALSE)
all.equal(t1, t3)  ## yay!
t2 <- test_sim(obs_disp = 100)

library(ggplot2)
(ggplot(t2, aes(date, value, colour = var))
  + geom_line()
  + geom_line(data=t1, lty=2, lwd=3, alpha=0.8)
)

library(rbenchmark)
b1 <- benchmark(test_sim(), test_sim(use_flex = FALSE),
          replications = 10)
## forecast_ensemble: not yet (needs calibrated model as input *or*
##  we need to change the interface)




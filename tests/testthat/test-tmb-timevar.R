Sys.setenv(R_TESTS="")

library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)

test_that('time variation works for several parameters and levels of continuity', {
  set_spec_version("0.1.1", "../../inst/tmb/")
  options(MP_use_state_rounding = FALSE)
  params <- read_params("ICU1.csv")
  test_fun = function(tv_dat) {
    mm = make_base_model(
      expand_params_S0(params, 1-1e-5),
      start_date = "2021-09-10",
      end_date = "2021-10-10",
      params_timevar = tv_dat,
      step_args = list(do_hazard = TRUE)
    )

    # change beta0, which is time-varying, so that
    # we can check that the parameter updates are
    # happening correctly in the C++ side
    test_pars = expand_params_S0(params, 1-1e-5)
    test_pars[1] = 1

    tmb_sim <- run_sim(
      params = test_pars,
      start_date = "2021-09-10",
      end_date = "2021-10-10",
      params_timevar = tv_dat,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      use_flex = TRUE,
      flexmodel = mm
    )

    r_sim <- run_sim(
      params = test_pars,
      start_date = "2021-09-10",
      end_date = "2021-10-10",
      params_timevar = tv_dat,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      use_flex = FALSE
    )

    compare_sims(r_sim, tmb_sim)
  }

  data.frame(
    Date = as.Date(c("2021-09-15", "2021-09-20", "2021-10-05", "2021-09-20", "2021-10-03")),
    Symbol = c("beta0", "beta0", "beta0", "Ca", "Ca"),
    Value = c(0.5, 0.1, 0.05, 0.1, 2.3),
    Type = c("rel_orig", "rel_orig", "rel_orig", "rel_orig", "rel_orig")
  ) %>% test_fun

  # FIXME: silence warning -- dropped switch times on final day
  data.frame(
    Date = ymd(20210910) + days(0:30),
    Symbol = "beta0",
    Value = seq(from = 0.1, to = 2, length.out = 31),
    Type = "rel_orig"
  ) %>% test_fun

  data.frame(
   Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
   Symbol = c("beta0", "beta0", "beta0"),
   Value = c(0.5, 0.1, 0.05),
   Type = c("rel_orig", "rel_orig", "rel_orig")
  ) %>% test_fun
})

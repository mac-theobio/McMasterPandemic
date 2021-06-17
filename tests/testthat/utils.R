## hack: stored object contains old copy of run_sim_breaks that uses
## old format for time_pars (Relative_value)
## FIXME: re-build and re-store
fix_stored <- function(x) {
  if (requireNamespace("anytime")) {
    x$forecast_args$sim_fun <- McMasterPandemic::run_sim_break
    x$forecast_args$start_date <- anytime::anydate(x$forecast_args$start_date)
    x$forecast_args$end_date <- anytime::anydate(x$forecast_args$end_date)
    if ("break_dates" %in% names(x$forecast_args$time_args)) {
      x$forecast_args$time_args$break_dates <-
        anytime::anydate(x$forecast_args$time_args$break_dates)
    }
  }
  return(x) ## can't convert dates, live with the warnings
}


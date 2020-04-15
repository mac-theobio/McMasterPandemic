% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrate.R
\name{forecast_sim}
\alias{forecast_sim}
\title{simulate based on a vector of parameters (including both time-varying change parameters, initial conditions, and other dynamical parameters), for fitting or forecasting}
\usage{
forecast_sim(
  p,
  opt_pars,
  base_params,
  start_date,
  end_date,
  break_dates,
  fixed_pars = NULL,
  sim_args = NULL,
  aggregate_args = NULL,
  return_val = c("aggsim", "vals_only")
)
}
\arguments{
\item{p}{vector of parameters}

\item{opt_pars}{starting parameters (and structure).  Parameters that are part of the \code{params_pansim} parameter vector can be specified within the \code{params} element (with prefixes if they are transformed); other parameters can include distributional parameters or time-varying parameters}

\item{base_params}{baseline parameters}

\item{start_date}{starting date for sims (far enough back to allow states to sort themselves out)}

\item{end_date}{ending date}

\item{break_dates}{specified breakpoints in beta0}

\item{fixed_pars}{parameters to fix}

\item{sim_args}{additional arguments to pass to \code{\link{run_sim}}}

\item{aggregate_args}{arguments passed to \code{\link{aggregate.pansim}}}

\item{return_val}{specify values to return (aggregated simulation, or just the values?)}
}
\description{
simulate based on a vector of parameters (including both time-varying change parameters, initial conditions, and other dynamical parameters), for fitting or forecasting
}
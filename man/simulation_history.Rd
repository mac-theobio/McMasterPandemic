% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexmodel_simulate.R
\name{simulation_history}
\alias{simulation_history}
\title{Simulation History}
\usage{
simulation_history(
  model,
  add_dates = TRUE,
  sim_params = NULL,
  include_initial_date = TRUE,
  obs_error = FALSE,
  condense = FALSE
)
}
\arguments{
\item{model}{a \code{\link{flexmodel}} object}

\item{add_dates}{should a column called \code{Date} be added?}

\item{sim_params}{vector to pass to a TMB AD function -- you can get
examples of these out of the \code{opt_par} element of a fitted
\code{\link{flexmodel}} or by calling \code{\link{tmb_params_trans}}}

\item{include_initial_date}{should the first row be the initial state
vector, or the first date after the initial?}

\item{obs_error}{should random observation error be added if it is
defined in the \code{model}? -- see \code{\link{update_error_dist}}}

\item{condense}{should a condensed set of simulation variables be
returned? -- see \code{\link{update_condense_map}}}
}
\description{
Simulation History
}
\seealso{
Other simulation: 
\code{\link{deprecated_simulation_functions}},
\code{\link{simulate_ensemble}()}
}
\concept{simulation}

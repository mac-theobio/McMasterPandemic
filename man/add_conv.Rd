% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmb.R
\name{add_conv}
\alias{add_conv}
\title{Add Variable by Convolution}
\usage{
add_conv(
  model,
  var_pattern,
  c_prop = "c_prop",
  c_delay_cv = "c_delay_cv",
  c_delay_mean = "c_delay_mean"
)
}
\arguments{
\item{model}{\code{\link{flexmodel}} object}

\item{var_pattern}{regular expression used to identify variables for
convolutions. See \code{\link{avail_for_conv}} for a function that will
return the names of all variables that are available for convolution.}

\item{c_prop}{name of the parameter in the \code{params} vector associated
with the proportion of individuals in the original state variable that are
represented in the convolved variable}

\item{c_delay_cv}{name of the parameter in the \code{params} vector
associated with the coefficient of variation of the gamma density}

\item{c_delay_mean}{name of the parameter in the \code{params} vector
associated with the mean of the gamma density}
}
\description{
Add a variable that will be available for comparisons with observed data
that is created by processing an existing variable in the simulation
history by convolving it with a gamma density.
}

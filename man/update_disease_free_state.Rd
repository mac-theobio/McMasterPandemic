% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmb.R
\name{update_disease_free_state}
\alias{update_disease_free_state}
\title{Update Disease-Free State}
\usage{
update_disease_free_state(model, state_pattern, param_pattern)
}
\arguments{
\item{model}{\code{\link{flexmodel}} object}

\item{state_pattern}{regular expression for identifying state variables
to be updated when constructing a disease-free state}

\item{param_pattern}{regular expression for identifying parameters
to use as disease-free state variables}
}
\value{
updated \code{\link{flexmodel}} object
}
\description{
Update Disease-Free State
}
\seealso{
Other flexmodel_definition_functions: 
\code{\link{add_factr}()},
\code{\link{add_linearized_outflow}()},
\code{\link{add_parallel_accumulators}()},
\code{\link{add_rate}()},
\code{\link{add_state_mappings}()},
\code{\link{add_state_param_sum}()},
\code{\link{initial_population}()},
\code{\link{mat_rate}()},
\code{\link{rep_rate}()},
\code{\link{update_linearized_params}()},
\code{\link{update_tmb_indices}()},
\code{\link{vec_rate}()}
}
\concept{flexmodel_definition_functions}

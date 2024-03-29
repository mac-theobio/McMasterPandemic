% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ageify.R
\name{expand_params_age}
\alias{expand_params_age}
\title{Expand parameter list to include age structure}
\usage{
expand_params_age(
  params,
  age_cat = mk_agecats(),
  beta0 = NULL,
  transmissibility = NULL,
  contact_rate_age = NULL,
  Nvec = NULL,
  pmat = NULL,
  balance_warning = TRUE
)
}
\arguments{
\item{params}{parameter list (e.g. read in with `read_params()`)}

\item{age_cat}{vector of age categories}

\item{beta0}{vector of age-specific beta0 values (if NULL, assume same beta0 for all age groups as already provided in params)}

\item{transmissibility}{proportion of contacts between S & I that lead to transmission}

\item{contact_rate_age}{average overall contact rate by age}

\item{Nvec}{population distribution (as counts); default is uniform}

\item{pmat}{contact matrix; default is uniform}

\item{balance_warning}{should a warning about the balance of contacts be provided?}
}
\description{
Expand parameter list to include age structure
}
\examples{
params <- read_params("PHAC_testify.csv")
params_age <- expand_params_age(params)
}
\seealso{
Other classic_macpan: 
\code{\link{add_d_log}()},
\code{\link{add_updated_vaxrate}()},
\code{\link{aggregate_agecats}()},
\code{\link{calibrate_comb}()},
\code{\link{calibrate}()},
\code{\link{check_age_cat_compatibility}()},
\code{\link{check_contact_rate_setting}()},
\code{\link{col_multiply}()},
\code{\link{condense_age}()},
\code{\link{condense_params_vax}()},
\code{\link{condense_state}()},
\code{\link{condense_vax}()},
\code{\link{dev_is_tikz}()},
\code{\link{do_step}()},
\code{\link{expand_params_desc_age}()},
\code{\link{expand_params_desc_variant}()},
\code{\link{expand_params_desc_vax}()},
\code{\link{expand_params_mistry}()},
\code{\link{expand_params_variant}()},
\code{\link{expand_params_vax}()},
\code{\link{expand_state_age}()},
\code{\link{expand_state_vax}()},
\code{\link{expand_stateval_testing}()},
\code{\link{fix_pars}()},
\code{\link{fix_stored}()},
\code{\link{forecast_ensemble}()},
\code{\link{forecast_sim}()},
\code{\link{getData}()},
\code{\link{get_GI_moments}()},
\code{\link{get_Gbar}()},
\code{\link{get_R0}()},
\code{\link{get_doses_per_day}()},
\code{\link{get_evec}()},
\code{\link{get_kernel_moments}()},
\code{\link{get_opt_pars}()},
\code{\link{get_r}()},
\code{\link{invlink_trans}()},
\code{\link{make_betavec}()},
\code{\link{make_beta}()},
\code{\link{make_jac}()},
\code{\link{make_ratemat}()},
\code{\link{make_state}()},
\code{\link{make_test_wtsvec}()},
\code{\link{make_vaxrate}()},
\code{\link{mk_Nvec}()},
\code{\link{mk_agecats}()},
\code{\link{mk_contact_rate_setting}()},
\code{\link{mk_mistry_Nvec}()},
\code{\link{mk_pmat}()},
\code{\link{mk_vaxcats}()},
\code{\link{mle_fun}()},
\code{\link{non_expanded_states}},
\code{\link{rExp}()},
\code{\link{read_params}()},
\code{\link{repair_names_age}()},
\code{\link{restore}()},
\code{\link{run_sim_ageify}()},
\code{\link{run_sim_break}()},
\code{\link{run_sim_loglin}()},
\code{\link{run_sim_mobility}()},
\code{\link{run_sim_range}()},
\code{\link{run_sim}()},
\code{\link{show_ratemat}()},
\code{\link{testify}()},
\code{\link{texify}()},
\code{\link{trans_state_vars}()},
\code{\link{update_contact_rate_setting}()},
\code{\link{update_foi}()},
\code{\link{update_params_mistry}()},
\code{\link{vis_model}()},
\code{\link{write_params}()}
}
\concept{classic_macpan}

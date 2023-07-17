#' Global options setup for using tmbstan with McMasterPandemic
#'
#' This is a hack: compiling the C++ code for tmbstan instead of using package-compiled objects
#'
#' @return NULL
#' @export
setup_stan = function(){
  cpp_dir = system.file('tmb', options()$MP_flex_spec_version, package = "McMasterPandemic")
  set_spec_version(options()$MP_flex_spec_version, cpp_dir, use_version_directories = FALSE)
}

#' Calibrate a flexmodel using STAN
#'
#' @param model a [McMasterPandemic::flexmodel_to_calibrate] object
#' @param chains number of MCMC chains
#'
#' @return a list with two elements
#'  - `model` the original model [McMasterPandemic::flexmodel_to_calibrate] object
#'  - `fit` the fit as a [rstan::stanfit] object
#' @export
calibrate_stan = function(
  model,
  model_to_calibrate,
  chains
){
  # extract tmb object
  model_tmb = tmb_fun(model_to_calibrate)

  fit = tmbstan::tmbstan(
    model_tmb,
    chains = chains
  )

  return(
    list(
      model = model, # need to carry this around for the forecast method because i can't get the simulation parameters to update for simuluation_history() using update_params() when the model object is a flexmodel_to_calibrate object... only seems to work with a flexmodel...
      # also simulation_history(as.flexmodel(flexmodel_to_calibrate)) crashes RStudio
      model_to_calibrate = model_to_calibrate,
      fit = fit
    )
  )
}

#' Tidy fit component output from [calibrate_stan()]
#'
#' @param model_fit the output of [calibrate_stan()]
#'
#' @return a list with two elements
#'  - `model` the original model [McMasterPandemic::flexmodel_to_calibrate] object
#'  - `fit` the fit as a [rstan::stanfit] object
tidy_fit_stan = function(
  model_fit
){
  McMasterPandemic::unpack(model_fit)

  # rename params
  names(fit) = c(
    names(tmb_params_trans(model_to_calibrate)),
    "log_posterior"
  )

  return(
    list(
      model = model,
      model_to_calibrate = model_to_calibrate,
      fit = fit
    )
  )
}

#' Produce a traceplot from a STAN fit
#'
#' @param model_fit the output of [calibrate_stan()]
#'
#' @return a traceplot of MCMC samples from STAN
#' @export
traceplot_stan = function(
    model_fit
){
  model_fit = tidy_fit_stan(model_fit)

  McMasterPandemic::unpack(model_fit)

  rstan::traceplot(fit, ncol = 1)
}

#' Forecast ensemble using STAN fit
#'
#' @param model_fit the output of [calibrate_stan()]
#' @param days_to_forecast the number of days to forecast in the future
#'
#' @return data.frame with individual realizations from ensemble
#' @export
#'
#' @examples
forecast_stan = function(
  model_fit,
  days_to_forecast,
  parallel = TRUE,
  n_cores = 2
){

  # unpack model and fit into local environment
  McMasterPandemic::unpack(model_fit)

  # collapse STAN samples into a single df
  samples = rstan::extract(fit)$params
  colnames(samples) = names(tmb_params_trans(model_to_calibrate))

  # simulate from STAN samples

  # register parallelization or sequential computation
  if(parallel){
    doParallel::registerDoParallel(cores = n_cores)
  } else{
    foreach::registerDoSEQ()
  }

  # extend simulation end date
  model = (model
    %>% McMasterPandemic::extend_end_date(
     days_to_extend = days_to_forecast
  ))

  # loop over samples, update model with sampled params, and simulate
  sims = foreach(i=c(1:nrow(samples)),
                 .packages=c('McMasterPandemic', 'magrittr')) %dopar% {
    params = samples[i,]

    # update model with samples and simulate
    # TODO: attach extension of end date
    # + time-varying forecast scenarios
    sim = McMasterPandemic::simulation_history(
     model
     # update model with sample params
     %>% McMasterPandemic::update_params(
       # back-transform any transformed values
       McMasterPandemic::invlink_trans(params)
     )
    )
  }

  # return list of raw iterations
  sims
}

#' Summarise forecast ensemble
#'
#' @param forecast output from [forecast_stan()]
#' @param var_order variable order in output (_e.g._ as output by [topological_sort()])
#' @param qvec named numeric vector giving quantiles to compute
#'
#' @return a data.frame with the ensemble summary
#' @export
summarise_forecast_stan = function(
    forecast,
    var_order,
    qvec = c(value = 0.5, lwr = 0.025, upr = 0.975)
){
  (forecast
    %>% summarise_trajectories(qvec = qvec)
    # remove rates
    %>% filter(!stringr::str_detect(var, "_to+_"))
    %>% mutate(var = factor(var, levels = var_order))
    %>% rename(date = Date)
    %>% arrange(date, var)
   )
}

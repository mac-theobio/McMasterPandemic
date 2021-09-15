#' Initialize Compartmental Model
#'
#' @param params a \code{param_pansim} object
#' @param state a \code{state_pansim} object
#' @param start_date simulation start date
#' @param end_date simulation end date
#' @param params_timevar data frame with scheduling for piece-wise
#' constant parameter variation (TODO: direct to other help pages)
#' @param do_hazard should hazard simulation steps be used?
#' (https://canmod.net/misc/flex_specs#v0.0.5) -- only used
#' if \code{spec_ver_gt('0.0.4')}
#' @return object representing a compartmental model
#' @export
init_model <- function(params, state,
                       start_date = NULL, end_date = NULL,
                       params_timevar = NULL,
                       do_hazard = TRUE, ...) {
  check_spec_ver_archived()
  model = list(
    state = state,
    params = params,
    ratemat = make_ratemat(state, params, sparse = TRUE),
    rates = list()
  )

  if(spec_ver_gt('0.0.1')) {
    model$parallel_accumulators = character()
  }

  if((!is.null(start_date)) & (!is.null(end_date))) {
    spec_check(introduced_version = '0.0.3', feature = 'Start and end dates')
    model$start_date = as.Date(start_date)
    model$end_date = as.Date(end_date)
    model$iters = as.integer(model$end_date - model$start_date)
  } else {
    if((!is.null(start_date)) | (!is.null(end_date))) {
      spec_check(introduced_version = '0.0.3', feature = 'Start and end dates')
      stop("\n\nIf you specify either a start or end date,\n",
           "you need to specify the other one as well.")
    }
    feature_check(introduced_version = '0.0.3', feature = 'Start and end dates')
  }

  if(spec_ver_gt('0.0.2')) {
    model$timevar = list(piece_wise = NULL)
  }
  if(!is.null(params_timevar)) {
    if(is.null(start_date) | is.null(end_date)) {
      stop("\n\nIf you specify a timevar table, you need to also\n",
           "specify a start and end date")
    }
    spec_check(introduced_version = '0.0.3',
               feature = 'Piece-wise contant time variation of parameters')


    if(spec_ver_eq('0.0.3')) {
      nbreaks = table(params_timevar$Symbol)
      pi_tv_par = find_vec_indices(names(nbreaks), model$params)
      names(pi_tv_par) = names(nbreaks)

      # want to be able to assume a particular structure
      # and ordering in downstream processing
      schedule = (
        params_timevar
        %>% mutate(Date = as.Date(Date))
        %>% mutate(Order = pi_tv_par[Symbol])
        %>% arrange(Order, Date)
      )

      o = order(pi_tv_par)
      pi_tv_par = pi_tv_par[o]
      nbreaks = nbreaks[o]

      # Initialize the values of the time-varying parameters
      # at the breakpoints
      #
      # this task is a bit like:
      # https://github.com/mac-theobio/McMasterPandemic/blob/478e520bf20279a2bd803da1301fc3cb15f03a34/R/sim_funs.R#L1050
      # the difference is that we are saving parameter values for
      # every breakpoint
      schedule$Init = NA
      new_param = TRUE
      ns = nrow(schedule)
      for(i in 1:ns) {
        if(new_param | schedule$Type[i] == 'rel_orig') {
          old_val = model$params[schedule$Symbol[i]]
        } else {
          old_val = schedule$Init[i-1]
        }
        schedule$Init[i] = old_val * schedule$Value[i]
        new_param = schedule$Symbol[i] != schedule$Symbol[min(ns, i+i)]
      }

      expanded_params = expand_time(model$params, schedule$Init, nbreaks, pi_tv_par)

      model$timevar$piece_wise = list(
        schedule = schedule,
        original_params = params,
        nbreaks = nbreaks,
        pi_tv_par = pi_tv_par)

      model$params = expanded_params
    } # v0.0.3

    if(spec_ver_gt("0.0.3")) {
      schedule = (
        params_timevar
        %>% mutate(Date = as.Date(Date))
        %>% mutate(breaks = as.integer(Date - model$start_date))
        %>% mutate(tv_spi = find_vec_indices(Symbol, c(state, params)))
        %>% arrange(breaks, tv_spi)
      )
      count_of_tv_at_breaks = c(table(schedule$breaks))

      # DRY: copied and modified from v0.0.3 above
      schedule$tv_val = NA
      new_param = TRUE
      ns = nrow(schedule)
      for(i in 1:ns) {
        if(new_param | schedule$Type[i] == 'rel_orig') {
          old_val = params[schedule$Symbol[i]]
        } else if(schedule$Type[i] == 'rel_prev') {
          old_val = schedule$tv_val[i-1]
        } else {
          stop("Unrecognized break-point type.\n",
               "Only rel_orig and rel_prev are allowed")
        }
        schedule$tv_val[i] = old_val * schedule$Value[i]
        new_param = schedule$Symbol[i] != schedule$Symbol[min(ns, i+i)]
      }

      attr(model$params, 'tv_param_indices') = setNames(find_vec_indices(
        unique(schedule$Symbol),
        params), unique(schedule$Symbol))


      model$timevar$piece_wise = list(
        schedule = schedule, # schedule includes tv_spi and tv_val as in spec
        breaks = as.integer(names(count_of_tv_at_breaks)),
        count_of_tv_at_breaks = unname(count_of_tv_at_breaks)
      )
    } # >v0.0.3
  }
  if(spec_ver_gt('0.0.4')) model$do_hazard = do_hazard

  # TODO: clarify index structure here once we converge
  model$tmb_indices = list()

  structure(model, class = "flexmodel")
}

#' Rate Structure
#'
#' Define how the rate of flow from one compartment to another
#' depends on the parameters and state variables.
#'
#' @param model compartmental model
#' @param from Name of state where flow is happening from
#' @param to Name of state where flow is happening to
#' @param formula Model formula defining dependence of the rate on
#' parameters and state variables
#' @return another compartmental model with an additional non-zero rate matrix
#' element specified
#' @export
add_rate = function(model, from, to, formula) {
  added_rate = (
    rate(from, to, formula, model$state, model$params, model$ratemat)
    %>% list
    %>% setNames(paste(from, to, sep = '_to_'))
  )
  model$rates = c(model$rates, added_rate)
  return(model)
}

#' @export
rate = function(from, to, formula, state, params, ratemat) {
  # TODO: test for formula structure
  # TODO: test that from and to are available in the state vector
  M = ratemat
  stopifnot(is.character(from),
            is.character(to),
            length(from) == 1L,
            length(to) == 1L,
            inherits(formula, 'formula'))

  # regex pattern for finding variables
  # (e.g. any parameter or state variable)
  # variable_regex looks like this '(beta0|Ca|...|zeta|S|E|Ia|...|V)'
  variable_regex = function(...) {
    character_class =
      (list(...)
       %>% lapply(names)
       %>% unlist
       %>% paste0(collapse = "|")
      )
    paste0('(', character_class, ')', sep = '')
  }
  get_variables = function(x) {
    r = regexpr(variable_regex(params, state), x)
    regmatches(x, r)
  }
  # FIXME: this only works because complements (1 - x) and
  # inverses (1 / x) are so similar in structure
  find_operators = function(x, operator) {
    grepl(
      paste0('\\( *1 *', operator, ' *', variable_regex(params, state),
             collapse = ''), x)
  }
  factor_table = function(x) {
    data.frame(
      var   = unlist(lapply(x, get_variables)),
      compl = unlist(lapply(x, find_operators, '-')),
      invrs = unlist(lapply(x, find_operators, '/')))
  }
  product_list = function(x) {
    x$factors = (x$formula
                 %>% as.character %>% getElement(2L)
                 %>% strsplit(split = '\\+') %>% getElement(1L)
                 %>% strsplit(split = '\\*')
                 %>% lapply(factor_table) %>% bind_rows(.id = 'prod_indx')
    )
    x$ratemat_indices =
      do.call(McMasterPandemic:::pfun, c(x[c('from', 'to')], list(mat = M)))
    x$factors$var_indx = find_vec_indices(x$factors$var, c(state, params))

    if(spec_ver_gt('0.0.1')) {
      x$state_dependent = any(x$factors$var_indx <= length(state))
    }
    if(spec_ver_gt('0.0.2')) {
      if(has_time_varying(params)) {
        x$factors$tv = x$factors$var %in%
          names(attributes(params)$tv_param_indices)
        x$time_varying = any(x$factors$tv)
      }
    }
    x
  }
  structure(
    product_list(list(from = from, to = to, formula = formula)),
    class = 'rate-struct')
}

find_vec_indices = function(x, vec) {
  (
    x
    %>% as.character
    %>% outer(names(vec), '==')
    %>% apply(1, which)
  )
}

#' Expand Parameters
#'
#' Should only be used for spec version 0.0.3
#'
#' @param x parameter vector
#' @param vals initial time-varying parameter values
#' @param nbreaks integer vector of time indices of the breakpoints
#' (in the case of piece-wise constant time variation, this can be
#' computed using \code{table(model$timevar$piece_wise$schedule$Symbol)})
#' @param pi_tv_par indices into the non-expanded/original parameter
#' vector for selecting time-varying parameters
#' @export
expand_time = function(x, vals, nbreaks, pi_tv_par){

  after = unname(pi_tv_par + c(0, nbreaks[-length(nbreaks)]))
  tv_param_indices = mapply(
    seq, from = after, length.out = nbreaks + 1, SIMPLIFY = FALSE)
  names(tv_param_indices) = names(nbreaks)

  names(vals) = paste0(
    rep(names(nbreaks), times = nbreaks),
    paste0('_t', sequence(nbreaks)))

  # insert time-dependent parameters
  start = 1
  end = 0
  for(i in seq_along(nbreaks)) {
    end = end + nbreaks[i]
    x = append(x, vals[start:end], after[i])
    start = end + 1
  }

  tv_param_indices = lapply(tv_param_indices, function(i) {
    setNames(i, names(x)[i])
  })

  structure(x, tv_param_indices = tv_param_indices)
}

#' @param x parameter vector
#' @export
has_time_varying = function(x) {
  spec_check(
    feature = "Time-varying parameters",
    introduced_version = '0.0.3')
  "tv_param_indices" %in% names(attributes(x))
}

#' @export
time_varying_rates = function(model) {
  model$rates[which_time_varying_rates(model)]
}

#' @export
which_time_varying_rates = function(model) {
  sd = sapply(model$rates, '[[', 'state_dependent')
  tv = sapply(model$rates, '[[', 'time_varying')
  which(sd | tv)
}

#' @export
state_dependent_rates = function(model) {
  model$rates[sapply(model$rates, '[[', 'state_dependent')]
}

#' Add Parallel Accumulators
#'
#' Add parallel accumulators to a compartmental model.
#'
#' @param model TODO
#' @param state_patterns regular expressions for identifying states as
#' parallel accumulators
#' @return another compartmental model with parallel accumulators specified
#' @export
add_parallel_accumulators = function(model, state_patterns) {
  if(spec_ver_gt('0.0.2')) {
    # TODO: check with experts if time-varying
    # parameters can be parallel accumulators?
    # currently assumed 'no' in the spec
    # https://canmod.net/misc/flex_specs#assumptions-0.0.3
    nms_tv_params = names(model$timevar$piece_wise$nbreaks)
    tv_accumulators = unlist(lapply(
      state_patterns, grep, nms_tv_params, perl = TRUE, value = TRUE))
    if(length(tv_accumulators) > 0) {
      stop('Time-varying parameters cannot be accumulators,\n',
           'but the following are:\n',
           paste(tv_accumulators, collapse = ', '))
    }
  }
  model$parallel_accumulators = parallel_accumulators(model, state_patterns)
  return(model)
}

#' @export
parallel_accumulators = function(model, state_patterns) {
  spec_check(introduced_version = '0.0.2', feature = "Parallel accumulators")
  (
    state_patterns
    %>% lapply(function(x) {grep(x, colnames(model$ratemat), value = TRUE)})
    %>% unlist
  )
}

#' Add TMB Indices
#'
#' Add, to a compartmental model, indices used to access appropriate values
#' during simulation and calibration using TMB
#'
#' @param model compartmental model
#' @param another compartmental model with indices for TMB
#' @export
add_tmb_indices = function(model) {
  model$tmb_indices = tmb_indices(model)
  return(model)
}

#' @export
ratemat_indices = function(rates, state_params) {
  sp = state_params
  ratemat_indices = sapply(rates, `[[`, 'ratemat_indices')
  spi = {lapply(rates, function(y) {y$factors$var_indx}) %>% unlist}
  count = sapply(rates, function(y) {
    nrow(y$factors)
  })
  modifier = lapply(unname(rates), '[[', 'factors') %>%
    bind_rows(.id = 'rate_indx') %>%
    mutate(add = as.logical(c(0, diff(as.numeric(prod_indx))))) %>%
    mutate(modifier = 4 * add + 2 * invrs + compl) %>%
    `$`("modifier")
  names(spi) = colnames(ratemat_indices) = names(count) = NULL
  indices = list(
    from = ratemat_indices[1,],
    to = ratemat_indices[2,],
    count = count,
    spi = spi,
    modifier = modifier
  )
  return(indices)
}

#' @export
piece_wise_indices = function(rates, timevar, start_date) {
  dd = lapply(unname(rates), '[[', 'factors') %>%
    bind_rows(.id = 'rate_indx')
  spi_tv_fac = with(dd, which(tv))
  nbreaks_fac = m$timevar$piece_wise$nbreaks[dd[spi_tv_fac,]$var]
  b = with(m$timevar$piece_wise$schedule, as.integer(Date - m$start_date))
  s = m$timevar$piece_wise$schedule$Symbol
  tbreaks = lapply(names(nbreaks_fac), function(v) {b[s == v]}) %>% unlist
}

#' @export
tmb_indices = function(model) {
  check_spec_ver_archived()

  sp = c(model$state, model$params)
  indices = list(make_ratemat_indices = ratemat_indices(model$rates, sp))

  if(spec_ver_gt('0.0.1')) {
    indices$par_accum_indices =
      which(colnames(model$ratemat) %in% model$parallel_accumulators)
  }
  if(spec_ver_eq('0.0.2')) {
    indices$update_ratemat_indices =
      ratemat_indices(state_dependent_rates(model), sp)
  }
  if(spec_ver_eq('0.0.3')) {
    indices$update_ratemat_indices =
      ratemat_indices(time_varying_rates(model), sp)
  }
  if(spec_ver_gt('0.0.3')) {
    indices$updateidx = which_time_varying_rates(model)
  }
  return(indices)
}

#' Make Objective Function with TMB
#'
#' Construct an objective function in TMB from a \code{flexmodel}
#' object. The behaviour of \code{tmb_fun} depends on \code{spec_version()}
#'
#' @param model object of class \code{flexmodel}
#' @param DLL argument to pass to the \code{DLL} argument of
#' \code{TMB::MakeADFun}. TODO: make this optional, and choose the main
#' package DLL once it exists
#' (https://github.com/mac-theobio/McMasterPandemic/issues/96)
#' @importFrom TMB MakeADFun
#' @useDynLib McMasterPandemic
#' @export
tmb_fun = function(model, DLL = "McMasterPandemic") {

  check_spec_ver_archived()

  # unpack model structure
  unpack(model)
  unpack(tmb_indices)
  unpack(make_ratemat_indices)
  if(spec_ver_gt('0.0.3')) unpack(timevar$piece_wise)


  # make ad functions
  if(spec_ver_eq('0.0.1')) {
    dd <- MakeADFun(data = list(state = c(state),
                                ratemat = ratemat,
                                from = from,
                                to = to,
                                count = count,
                                spi = spi,
                                modifier = modifier),
                    parameters = list(params=c(params)),
                    DLL=DLL)
  } else if(spec_ver_eq('0.0.2')) {
    up = update_ratemat_indices

    # TODO: spec problem -- numIters should have been kept in model
    numIters = 3

    dd <- MakeADFun(data = list(state = c(state),
                                ratemat = ratemat,
                                from = from,
                                to = to,
                                count = count,
                                spi = spi,
                                modifier = modifier,
                                update_from = up$from,
                                update_to = up$to,
                                update_count = up$count,
                                update_spi = up$spi,
                                update_modifier = up$modifier,
                                par_accum_indices = par_accum_indices,
                                numIterations = numIters),
                    parameters = list(params=c(params)),
                    DLL=DLL)
  } else if(spec_ver_eq('0.0.4')) {

    dd <- MakeADFun(data = list(state = c(state),
                                ratemat = ratemat,
                                from = from,
                                to = to,
                                count = count,
                                spi = spi,
                                modifier = modifier,
                                updateidx = c(updateidx),
                                breaks = breaks,
                                count_of_tv_at_breaks = count_of_tv_at_breaks,
                                tv_spi = schedule$tv_spi,
                                tv_val = schedule$tv_val,
                                par_accum_indices = par_accum_indices,
                                numIterations = iters),
                    parameters = list(params=c(params)),
                    DLL=DLL)

  } else if(spec_ver_eq('0.0.5')) {

    dd <- MakeADFun(data = list(state = c(state),
                                ratemat = ratemat,
                                from = from,
                                to = to,
                                count = count,
                                spi = spi,
                                modifier = modifier,
                                updateidx = c(updateidx),
                                breaks = breaks,
                                count_of_tv_at_breaks = count_of_tv_at_breaks,
                                tv_spi = schedule$tv_spi,
                                tv_val = schedule$tv_val,
                                par_accum_indices = par_accum_indices,
                                do_hazard = do_hazard,
                                numIterations = iters),
                    parameters = list(params=c(params)),
                    DLL=DLL)
  } else{
    stop('This feature is not supported by your installation of MacPan')
  }
  return(dd)
}

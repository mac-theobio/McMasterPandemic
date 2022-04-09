##' Initialize Flexible Compartmental Model
##'
##' https://canmod.net/misc/flex_specs
##'
##' @param params a \code{param_pansim} object
##' @param state a \code{state_pansim} object
##' @param start_date simulation start date
##' @param end_date simulation end date
##' @param params_timevar data frame with scheduling for piece-wise
##' constant parameter variation (TODO: direct to other help pages)
##' @param do_hazard should hazard simulation steps be used?
##' (https://canmod.net/misc/flex_specs#v0.0.5) -- only used
##' if \code{spec_ver_gt('0.0.4')}
##' @param do_hazard_lin like \code{do_hazard} but for the
##' linearized model that is used to construct the initial state
##' variable -- only used when \code{do_make_state == TRUE}
##' @param do_approx_hazard approximate the hazard transformation
##' by a smooth function (experimental)
##' @param do_approx_hazard like \code{do_approx_hazard} but for
##' the linearized model that is used to construct the initial
##' state (experimental)
##' @param do_make_state should state be remade on the c++ size?
##' (https://canmod.net/misc/flex_specs#v0.1.1) -- only used
##' if \code{spec_ver_gt('0.1.0')}
##' @param max_iters_eig_pow_meth maximum number of iterations
##' to use in computing the eigenvector for initial state
##' construction
##' @param tol_eig_pow_meth tolerance for determining convergence
##' of the power method used in initial state construction
##' @param data optional observed data frame in long format to
##' compare with simulated trajectories. must have the following
##' columns: \code{date}, \code{var}, \code{value}.
##' @family flexmodels
##' @return flexmodel object representing a compartmental model
##' @importFrom lubridate Date
##' @export
flexmodel <- function(params, state = NULL,
                       start_date = NULL, end_date = NULL,
                       params_timevar = NULL,
                       do_hazard = getOption("MP_default_do_hazard"),
                       do_hazard_lin = FALSE,
                       do_approx_hazard = FALSE,
                       do_approx_hazard_lin = TRUE,
                       do_make_state = getOption("MP_default_do_make_state"),
                       max_iters_eig_pow_meth = 8000,
                       tol_eig_pow_meth = 1e-6,
                       haz_eps = 1e-6,
                       data = NULL,
                       ...) {
    check_spec_ver_archived()
    name_regex = "^" %+% getOption("MP_name_search_regex") %+% "$"
    if(!all(grepl(name_regex, c(names(params), names(state))))) {
        stop("only syntactically valid r names can be used ",
             "for state variables and parameters")
    }

    if (is.null(state)) {
      if (!inherits(params, "params_pansim")) {
        stop("an initial state vector is required, because\n",
             "params is not of class params_pansim")
      }
      # inefficient! should just directly make a zero'd state vector.
      # trying to be more efficient:
      #  - tried setting use_eigvec = FALSE, but this failed for some reason (bug??)
      #  - for now we can do this ugly thing of turning down the number of power
      #    method steps and then restoring
      n_steps_default = getOption("MP_rexp_steps_default")
      options(MP_rexp_steps_default = 1)
      state = make_state(params = params)
      options(MP_rexp_steps_default = n_steps_default)
      state[] = 0
    } else if (is.character(state)) {
      state = setNames(rep(0, length(state)), state)
    }

    if(inherits(state, "state_pansim") & inherits(params, "params_pansim")) {
      ratemat = make_ratemat(state, params, sparse = TRUE)
    } else {
      ratemat = matrix(
        0,
        nrow = length(state), ncol = length(state),
        dimnames = list(from = names(state), to = names(state))
      )
    }

    # TODO: keep an eye on this -- i think that the
    # flex framework should _not_ use parameter
    # lists and rather stick to numeric vectors, but
    # not totally sure
    # -- also should use get_attr and put_attr from utils.R
    pattr = attributes(params)
    params = setNames(unlist(params), names(params))  # FIXME: will silently fail for nested lists
    attributes(params) = c(attributes(params), pattr)

    model <- list(
        state = state,
        params = params,
        ratemat = ratemat,
        rates = list(),
        name_regex = name_regex,
        spec_ver = as.character(spec_version())
    )

    if (spec_ver_btwn("0.0.1", "0.1.1")) {
        # parallel accumulator declaration now handled by outflow
        model$parallel_accumulators <- character(0L)
    }

    if ((!is.null(start_date)) & (!is.null(end_date))) {
        spec_check(introduced_version = "0.0.3",
                   feature = "Start and end dates")
        model$start_date <- as.Date(start_date)
        model$end_date <- as.Date(end_date)
        model$iters <- (model$end_date
            %>% difftime(model$start_date, units = 'days')
            %>% as.integer
        )
        if (model$iters < 0) {
          stop("start_date must be less than or equal to end_date")
        }
    } else {
        if ((!is.null(start_date)) | (!is.null(end_date))) {
            spec_check(introduced_version = "0.0.3",
                       feature = "Start and end dates")
            stop(
                "\n\nIf you specify either a start or end date,\n",
                "you need to specify the other one as well."
            )
        }
        feature_check(introduced_version = "0.0.3",
                      feature = "Start and end dates")
    }

    if (spec_ver_gt("0.0.2")) {
        model$timevar <- list(piece_wise = NULL)
    }
    if (!is.null(params_timevar)) {
        if (is.null(start_date) | is.null(end_date)) {
            stop(
                "\n\nIf you specify a timevar table, you need to also\n",
                "specify a start and end date"
            )
        }
        spec_check(
            introduced_version = "0.0.3",
            feature = "Piece-wise contant time variation of parameters"
        )

        if (spec_ver_gt("0.0.3")) {
          model = update_piece_wise(model, params_timevar)
        } ## >v0.0.3
    } else {
        model = initialize_piece_wise(model)
    }

    if (spec_ver_gt("0.0.4")) model$do_hazard <- do_hazard

    if (spec_ver_gt("0.0.6")) model$sums = list()
    model$sum_vector = numeric(0L)

    if (spec_ver_gt("0.1.0")) {
        model$do_hazard_lin <- do_hazard_lin
        model$do_approx_hazard = do_approx_hazard
        model$do_approx_hazard_lin = do_approx_hazard_lin
        if(max_iters_eig_pow_meth < 100) {
          warning("maximum number of iterations for the power method must be at least 100 -- setting max_iters_eig_pow_meth = 100")
          max_iters_eig_pow_meth = 100L
        }
        model$do_make_state = do_make_state
        model$max_iters_eig_pow_meth = max_iters_eig_pow_meth
        model$tol_eig_pow_meth = tol_eig_pow_meth
        model$haz_eps = haz_eps
        model$disease_free = list()
        model$linearized_params = list()
        model$outflow = list()
        model$linearized_outflow = list()

        model$initialization_mapping = init_initialization_mapping
        model$initial_population = init_initial_population

        which_step_zero_tv = which(
          model$timevar$piece_wise$breaks == 0L)
        model$timevar$piece_wise$step_zero_tv_idx =
          model$timevar$piece_wise$schedule$tv_spi[which_step_zero_tv]
        model$timevar$piece_wise$step_zero_tv_vals =
          model$timevar$piece_wise$schedule$tv_val[which_step_zero_tv]
        model$timevar$piece_wise$step_zero_tv_count = rep(1, length(which_step_zero_tv))
    }

    if (spec_ver_gt("0.1.1")) {
      model$factrs = list()
      model$sim_report_exprs = list()
    }
    model$factr_vector = init_factr_vector

    # condensation -- spec_ver_gt("0.1.2)
    model$condensation = init_condensation

    if (spec_ver_gt("0.1.2") & !is.null(data)) {
      stopifnot(isTRUE(all.equal(c(names(data)), c("date", "var", "value"))))
      #allvars = model$condensation_map[final_sim_report_names(model)]
      obsvars = unique(data$var)
      #stopifnot(all(obsvars %in% allvars))
      model$observed$data = data

      # in the future, the user should be able to provide this
      # loss_params data themselves. would open up the
      # possibility to add error distributions other than
      # the negative binomial, including distributions with
      # more than one parameter (in addition to the location
      # parameter that is determined by the simulations).
      # right now the c++ side assumes only negative binomial.
      model$observed$loss_params = data.frame(
        Parameter = "nb_disp", # only choice: dispersion
        Distribution = "nb",   # only choice: negative binomial
        Variable = obsvars
      )
      model$params = expand_params_nb_disp(model$params, obsvars)
    } else {
      model$observed = init_observed
    }

    model$opt_params = list()
    model$opt_tv_params = list()
    model$condensation_map = init_condensation_map

    if (FALSE & spec_ver_gt('0.1.2')) {
      model$opt_params = initialize_opt_params(model)
      model$opt_tv_params = initialize_opt_tv_params(model)
    }

    model$tmb_indices <- init_tmb_indices

    structure(model, class = "flexmodel")
}

#' `init_model` is deprecated and identical to `flexmodel`.
#'
#' @rdname flexmodel
#' @export
init_model = function(...) {
  warning("\ninit_model is deprecated.\nplease use flexmodel.\n",
          "the only difference between the two functions is the name.")
  flexmodel(...)
}

# rate and associated functions ---------------------
#
#   add_rate
#   rep_rate
#   vec_rate
#   mat_rate

##' Define Rate for Single Element of Rate Matrix
##'
##' @param from from state
##' @param to to state
##' @param formula one-sided formula defining the rate with reference
##' to the parameters and state variables
##' @param state state_pansim object
##' @param params param_pansim object
##' @param sums vector of sums of state variables and parameters
##' @param ratemat rate matrix
##' @importFrom dplyr bind_rows
##' @family flexmodels
##' @export
rate <- function(from, to, formula, state, params, sums, factrs, ratemat) {
    ## TODO: test for formula structure
    M <- ratemat
    stopifnot(
        is.character(from),
        is.character(to),
        length(from) == 1L,
        length(to) == 1L,
        (  inherits(formula, "formula")
         | is.character(formula)
         | is(formula, 'struc')
        )
    )
    if (!from %in% names(state)) {
      stop("state variable ", from, " used but not found in the model")
    }
    if (!to %in% names(state)) {
      stop("state variable ", to, " used but not found in the model")
    }

    product_list <- function(x) {
        x$factors <- (x$formula
            %>% parse_formula
            %>% lapply(factor_table) %>% bind_rows(.id = "prod_indx")
        )
        x$ratemat_indices <- do.call(
            pfun,
            c(x[c("from", "to")], list(mat = M)))
        if(nrow(x$ratemat_indices) > 1L) {
            stop('More than one element of the rate matrix is being ',
                 'referred to.\nTry using rep_rate instead of rate.')
        }
        x$factors$var_indx <- find_vec_indices(
            x$factors$var,
            c(state, params, sums, factrs))

        missing_vars = x$factors$var[sapply(x$factors$var_indx, length) == 0L]
        if(length(missing_vars) > 0L) {
            stop("The following variables were used to define the model,\n",
                 "but they could not be found in the state, parameter or sum ",
                 "vectors:\n",
                 paste0(missing_vars, collapse = "\n"))
        }

        if (spec_ver_gt("0.0.1")) {
            x$state_dependent <- any(x$factors$var_indx <= length(state))
        }
        if (spec_ver_gt("0.0.2")) {
            if (has_time_varying(params)) {
                x$factors$tv <- x$factors$var %in%
                    names(attributes(params)$tv_param_indices)
                x$time_varying <- any(x$factors$tv)
            } else {
                x$time_varying = FALSE
            }
        }
        if (spec_ver_gt('0.0.6')) {
            x$sum_dependent =
                any(x$factors$var_indx > (length(state) + length(params)))
        } else {
            x$sum_dependent = FALSE
        }
        x
    }

    structure(
        product_list(list(from = from, to = to, formula = formula)),
        class = "rate-struct"
    )
}

##' Rate Structure
##'
##' Define how the rate of flow from one compartment to another
##' depends on the parameters and state variables.
##'
##' @param model compartmental model
##' @param from Name of state where flow is happening from
##' @param to Name of state where flow is happening to
##' @param formula Model formula defining dependence of the rate on
##' parameters and state variables
##' @return another compartmental model with an additional non-zero rate matrix
##' element specified
##' @family flexmodels
##' @export
add_rate <- function(model, from, to, formula) {
    unpack(model)
    added_rate <- (from
      %>% rate(to, formula, state, params, sum_vector, factr_vector, ratemat)
      %>% list
      %>% setNames(paste(from, to, sep = "_to_"))
    )
    model$rates <- c(model$rates, added_rate)
    return(model)
}

#' Repeat a Rate for Several Rate Matrix Elements
#'
#' @param model flexmodel
#' @param formula formula or length-1 character vector
#' @family flexmodels
#' @export
rep_rate = function(model, from, to, formula,
                    mapping = c("pairwise", "blockwise")) {

    map_fun = switch(
        match.arg(mapping),
        pairwise = pwise,
        blockwise = block)

    stopifnot(inherits(model, "flexmodel"))

    unpack(model)
    #check_from_to(from, to, names(state))
    indices = map_fun(from, to, model$ratemat)

    if(!inherits(indices, "matrix")) {
        stop("indices must be a matrix")
    } else if(!(ncol(indices) == 2L)) {
        stop("indices must be a two-column matrix")
    } else if(!all(colnames(indices) == c("from_pos", "to_pos"))) {
        stop("indices must be a matrix with columns from_pos and to_pos")
    }

    if(!inherits(formula, "formula")) {
        if(!inherits(formula, "character")) {
            stop("formula must be either a formula or character object")
        } else if(length(formula) != 1L){
            stop("character formulas must be of length 1")
        }
    }

    from = rownames(ratemat)[indices[,'from_pos']]
    to = colnames(ratemat)[indices[,'to_pos']]

    lst = mapply(rate, from, to,
        MoreArgs = nlist(formula, state, params, sums, factrs = factr_vector, ratemat),
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    nms = mapply(paste, from, to, MoreArgs = list(sep = "_to_"))
    model$rates <- c(rates, setNames(lst, nms))

    return(model)
}

#' Specify Vector of Rates
#'
#' @family flexmodels
#' @export
vec_rate = function(model, from, to, formula,
                    mapping = c("pairwise", "blockwise")) {

    unpack(model)
    #check_from_to(from, to, names(state))
    map_fun = switch(
        match.arg(mapping),
        pairwise = pwise,
        blockwise = block)

    indices = map_fun(from, to, ratemat)

    from = rownames(ratemat)[indices[,'from_pos']]
    to = colnames(ratemat)[indices[,'to_pos']]

    lst = mapply(rate, from, to, as.character(formula),
                 MoreArgs = nlist(state, params, sums, factrs = factr_vector, ratemat),
                 SIMPLIFY = FALSE, USE.NAMES = FALSE)
    nms = mapply(paste, from, to, MoreArgs = list(sep = "_to_"))
    model$rates <- c(rates, setNames(lst, nms))

    return(model)
}

#' Specify Matrix of Rates
#'
#' Not implemented
#'
#' @family flexmodels
#' @export
mat_rate = function() {
    stop("\nrate specification with matrices is ",
         "coming sometime in the future ... maybe\n",
         "in the meantime you can specify vector-valued rates with vec_rate\n",
         "see this document for potentially more information on priorities:\n",
         options("MP_flex_spec_doc_site")[[1]])
}



# factr and associated functions ----------------------

#' @export
factr <- function(factr_nm, formula, state, params, sums, factrs, ratemat) {
  ## TODO: test for formula structure
  stopifnot(
    (  inherits(formula, "formula")
       | is.character(formula)
       | is(formula, 'struc')
    )
  )

  product_list <- function(x) {
    factor_table = factor_table
    find_vec_indices = find_vec_indices
    spec_check(
      introduced_version = "0.1.2",
      feature = "common factors (i.e. factr)"
    )
    x$factors <- (x$formula
                  %>% parse_formula
                  %>% lapply(factor_table) %>% bind_rows(.id = "prod_indx")
    )

    x$factors$var_indx <- find_vec_indices(
      x$factors$var,
      c(state, params, sums, factrs))

    missing_vars = x$factors$var[sapply(x$factors$var_indx, length) == 0L]
    if(length(missing_vars) > 0L) {
      stop("The following variables were used to define the model,\n",
           "but they could not be found in the state, parameter or sum ",
           "vectors:\n",
           paste0(missing_vars, collapse = "\n"))
    }

    x$state_dependent <- any(x$factors$var_indx <= length(state))
    if (has_time_varying(params)) {
      x$factors$tv <- x$factors$var %in%
        names(attributes(params)$tv_param_indices)
      x$time_varying <- any(x$factors$tv)
    } else {
      x$time_varying = FALSE
    }

    x$sum_dependent = any(x$factors$var_indx > (length(state) + length(params)))
    x
  }

  structure(
    product_list(list(factr_nm = factr_nm, formula = formula)),
    class = "factr-struct"
  )
}

#' Intermediate Factors
#'
#' Save intermediate computations so that they can be used
#' as factors in rate expressions and returned as part
#' of the simulation history. The definition of a factor
#' is given here: \url{https://canmod.net/misc/flex_specs#v0.1.0}.
#'
#' @param model \code{\link{flexmodel}} object
#' @param factr_nm name of the new intermediate factor
#' @param formula formula (or string) that follows
#' this spec, \url{https://canmod.net/misc/flex_specs#v0.1.0},
#' or 1-by-1 \code{\link{struc}} object describing the
#' intermediate factor
#'
#' @seealso See \code{\link{vec_factr}} to add more than
#' one intermediate factor at the same time.
#'
#' @return \code{\link{flexmodel}} object
#' @family flexmodels
#' @export
add_factr <- function(model, factr_nm, formula) {
  unpack(model)

  added_factr <- (factr_nm
                 %>% factr(formula,
                           state, params,
                           sum_vector, factr_vector)
                 %>% list
                 %>% setNames(factr_nm)
  )

  added_factr[[factr_nm]]$initial_value = eval_formulas(list(formula), get_var_list(model))
  model$factrs <- c(model$factrs, added_factr)

  model$factr_vector = get_factr_initial_value(model) %>% unlist
  return(model)
}

#' Vectors of Intermediate Factors
#'
#' @param model \code{\link{flexmodel}} object
#' @param factr_nms vector of names of the new intermediate factors
#' @param formula \code{\link{struc}} object describing the
#' vector of intermediate factors
#'
#' @seealso See \code{\link{add_factr}} to add a single
#' scalar-valued intermediate factor and \code{\link{vec}}
#' to create a vector-valued \code{\link{struc}} object.
#'
#' @return \code{\link{flexmodel}} object
#'
#' @export
vec_factr = function(model, factr_nms, formula) {

  unpack(model)
  lst = mapply(factr, factr_nms, as.character(formula),
               MoreArgs = nlist(state, params, sums, factrs = factr_vector),
               SIMPLIFY = FALSE, USE.NAMES = FALSE)

  for(i in seq_along(lst)) {
    lst[[i]]$initial_value = eval_formulas(list(lst[[i]]$formula), get_var_list(model))
  }
  model$factrs <- c(factrs, setNames(lst, factr_nms))

  model$factr_vector = get_factr_initial_value(model) %>% unlist
  return(model)
}

# sums of state variables and parameters ---------------------

state_param_sum = function(sum_name, summands, state, params) {
    spec_check('0.1.0', 'sums of state variables and parameters')
    if (!is.character(sum_name)) stop("sum_name must be character-valued")
    if (length(sum_name) != 1L) {
      stop("can only specify one sum at a time, ",
           "but sum_name was a vector with ", length(sum_name),
           " elements.")
    }
    sp = c(state, params)
    if (sum_name %in% names(sp)) {
      stop("sums cannot be named after state variables ",
           "or parameters.")
    }
    summands = (summands
        %>% lapply(grep, names(sp), value = TRUE)
        %>% unlist
    )
    if (length(summands) == 0L) {
      stop("regular expressions did not match any ",
           "state variables or parameters to sum ",
           "and store as ", sum_name)
    }
    ii = find_vec_indices(summands, sp)
    val = sum(sp[ii])
    list(
        summands = summands,
        sum_indices = ii,
        initial_value = val)
}

#' Sums of States and Parameters
#'
#' Save intermediate sums of states and parameters so that
#' they can be used as factors in rate expressions and
#' retured as part of the simulation history.
#'
#' @param model \code{\link{flexmodel}}
#' @param sum_name name of sum of state variables and parameters
#' @param summands character vector of regular expressions for identifying
#' state variables and parameters to sum together
#'
#' @return \code{\link{flexmodel}}
#'
#' @family flexmodels
#' @export
add_state_param_sum = function(model, sum_name, summands) {
    if (length(model$factrs) != 0L) {
      stop("cannot add any more state-param sums after intermediate factrs have been added")
    }
    model$sums[[sum_name]] = state_param_sum(
        sum_name, summands, model$state, model$params)

    # assumes that order of sums doesn't change!
    model$sum_vector = get_sum_initial_value(model) %>% unlist
    model
}


# condensation -----------------------------------

#' @export
sim_report_expr = function(expr_nm, formula, init_sim_report_nms) {
  stopifnot(
    (  inherits(formula, "formula")
       | is.character(formula)
       | is(formula, 'struc')
    )
  )

  product_list <- function(x) {
    factor_table = factor_table
    find_vec_indices = find_vec_indices
    spec_check(
      introduced_version = "0.1.2",
      feature = "simulation reports"
    )
    x$factors <- (x$formula
                  %>% parse_formula
                  %>% lapply(factor_table) %>% bind_rows(.id = "prod_indx")
    )

    x$factors$var_indx <- find_vec_indices(
      x$factors$var,
      init_sim_report_nms)

    # TODO -- add this check copied from factr back for sim_report_expr
    # missing_vars = x$factors$var[sapply(x$factors$var_indx, length) == 0L]
    #if(length(missing_vars) > 0L) {
    #  stop("The following variables were used to define the model,\n",
    #       "but they could not be found in the state, parameter or sum ",
    #       "vectors:\n",
    #       paste0(missing_vars, collapse = "\n"))
    #}

    # This stuff about dependence shouldn't matter for sim_report_expr
    # because everything changes with time -- that's kind of the point
    #
    # x$state_dependent <- any(x$factors$var_indx <= length(state))
    # if (has_time_varying(params)) {
    #   x$factors$tv <- x$factors$var %in%
    #     names(attributes(params)$tv_param_indices)
    #   x$time_varying <- any(x$factors$tv)
    # } else {
    #   x$time_varying = FALSE
    # }
    # x$sum_dependent = any(x$factors$var_indx > (length(state) + length(params)))

    x
  }

  structure(
    product_list(list(expr_nm = expr_nm, formula = formula)),
    class = "sim-report-expr-struct"
  )

}

#' @export
add_sim_report_expr = function(model, expr_nm, formula) {
  unpack(model)

  added_expr <- (expr_nm
    %>% sim_report_expr(
      formula,
      initial_sim_report_names(model)
    )
    %>% list
    %>% setNames(expr_nm)
  )

  # I don't think we need any initialization of these values
  # TODO -- verify this
  # added_expr[[expr_nm]]$initial_value = eval_formulas(list(formula), get_var_list(model))
  model$sim_report_exprs = c(model$sim_report_exprs, added_expr)
  # model$factr_vector = get_factr_initial_value(model) %>% unlist

  return(model)
}

#' @export
add_lag_diff = function(
  model, var_pattern,
  delay_n = 1) {
  stopifnot(is_len1_int(delay_n))
  stopifnot(is_len1_char(var_pattern))
  var_matches = grep(var_pattern
    , intermediate_sim_report_names(model)
    , perl = TRUE
    , value = TRUE
  )
  output_names = "lag" %_% delay_n %_% "diff" %_% var_matches
  added_lag_diff = nlist(
    var_pattern,
    delay_n,
    output_names
  )
  model$lag_diff = c(
    model$lag_diff,
    list(added_lag_diff)
  )
  model
}

#' @export
add_conv = function(
  model, var_pattern,
  c_prop = "c_prop",
  c_delay_cv = "c_delay_cv",
  c_delay_mean = "c_delay_mean") {

  stopifnot(is_len1_char(var_pattern))

  var_matches = grep(var_pattern
    , intermediate_sim_report_names(model)
    , perl = TRUE
    , value = TRUE
  )
  output_names = "conv" %_% var_matches

  added_conv = list(
    var_pattern = var_pattern,
    conv_pars = nlist(c_prop, c_delay_cv, c_delay_mean),
    output_names = output_names
  )
  model$conv = c(
    model$conv,
    list(added_conv)
  )
  model
}

#' Specify and Name Variables for Condensation
#'
#' The condensed simulation is a (possibly renamed) subset
#' of the variables in the simulation history. This function
#' allows one to define this subset and how it is renamed,
#' by creating a map with the following structure:
#' \code{c(orig_var_i = "cond_var_1", ..., orig_var_j = "cond_var_n")}.
#'
#' @param model flexmodel
#' @param map named vector with names that are a subset of
#' the variables in the simulation model. if \code{NULL} the
#' identity map is used that makes all simulation history
#' variables available without name changes.
#' (\code{final_sim_report_names(model)}) and values that
#' give the names of the variables in the condensed data set
#' @export
update_condense_map = function(model, map = NULL) {
  allvars = final_sim_report_names(model)
  if (is.null(map)) map = setNames(allvars, allvars)
  stopifnot(all(names(map) %in% allvars))
  model$condensation_map = map
  model
}


# parallel accumulators (deprecated -- use outflow instead) ---------------------

##' Add Parallel Accumulators
##'
##' Add parallel accumulators to a compartmental model.
##'
##' @param model TODO
##' @param state_patterns regular expressions for identifying states as
##' parallel accumulators
##' @return another compartmental model with parallel accumulators specified
##' @family flexmodels
##' @export
add_parallel_accumulators <- function(model, state_patterns) {
    model$parallel_accumulators <- parallel_accumulators(model, state_patterns)
    return(model)
}

##' @family flexmodels
##' @export
parallel_accumulators <- function(model, state_patterns) {
    spec_check(introduced_version = "0.0.2", feature = "Parallel accumulators")
    if(spec_ver_gt('0.1.0')) stop('Parallel accumulators are now handled through outflow')
    (state_patterns
        %>% lapply(function(x) {
            grep(x, colnames(model$ratemat), value = TRUE)
        })
    )
}


##' @rdname add_outflow
##' @export
add_linearized_outflow = function(model, from, to) {
    spec_check(
      introduced_version = '0.1.1',
      feature = 'Flexible restriction of outflows in the linearized model'
    )
    model$linearized_outflow = append(
        model$linearized_outflow,
        list(outflow(model, from, to)))
    return(model)
}

##' Add Outflows
##'
##' Add outflows corresponding to inflows specified by
##' \code{add_rate}, \code{rep_rate}, and \code{vec_rate}.
##'
##' By default, this function will set outflows for all
##' corresponding inflows. To define outflows \code{from}
##' specific states \code{to} other specific states, pass
##' regular expressions to the \code{from} and \code{to}
##' arguments to identify these states.
##'
##' \code{add_linearized_outflow} is used to specify
##' outflows in linearized models
##' (\url{https://canmod.net/misc/flex_specs#v0.1.1}).
##'
##' @param model \code{\link{flexmodel}} object
##' @param from string giving a regular expression for identifying
##' states from which individuals are flowing
##' @param to string giving a regular expression for identifying
##' states to which individuals are flowing
##'
##' @return model \code{\link{flexmodel}} object
##'
##' @family flexmodels
##' @export
add_outflow = function(
  model,
  from = '.+',
  to = '.+') {

  model$outflow = append(
    model$outflow,
    list(outflow(model, from, to)))
  return(model)
}

##' @family flexmodels
##' @export
outflow = function(
  model,
  from = '.+',
  to = '.+') {

  spec_check(
    introduced_version = '0.1.1',
    feature = 'Flexible restriction of outflows'
  )
  nlist(from, to)
}

##' @family flexmodels
##' @export
update_linearized_params = function(model, param_pattern, value) {
    model$linearized_params = c(
        model$linearized_params,
        list(linearized_params(model, param_pattern, value)))
    return(model)
}

##' @family flexmodels
##' @export
linearized_params = function(model, param_pattern, value) {
    spec_check(introduced_version = "0.1.1",
               feature = "Disease free parameter updates")
    params_to_update = grep(param_pattern,
                            names(model$params),
                            value = TRUE, perl = TRUE)
    if(length(params_to_update) == 0L)
        stop("param_pattern does not match any parameters")
    update_value = value
    nlist(params_to_update, update_value)
}

##' @family flexmodels
##' @export
update_disease_free_state = function(model, state_pattern, param_pattern) {
    model$disease_free = c(
        model$disease_free,
        list(disease_free_state(model, state_pattern, param_pattern)))
    return(model)
}

##' @family flexmodels
##' @export
disease_free_state = function(model, state_pattern, param_pattern) {
    spec_check(introduced_version = "0.1.1",
               feature = "Disease free state updates")
    states_to_update = grep(state_pattern,
                            names(model$state),
                            value = TRUE, perl = TRUE)
    params_to_use = grep(param_pattern,
                         names(model$params),
                         value = TRUE, perl = TRUE)
    if(length(params_to_use) != 1L) {
        if(length(states_to_update) != length(params_to_use)) {
            stop("incompatible mapping of parameters to states")
        }
    }
    nlist(states_to_update, params_to_use)
}

##' Initial Population
##'
##' \code{infected} is the multiplier of the normalized infected components of the
##' eigenvector. \code{total - infected} is multiplied by
##' \code{1/length(initial_susceptible)}
##'
##' In the future we might want the ability to specify a distribution
##' across susceptible compartments -- currently a uniform distribution
##' is assumed..
##'
##' @param total name of a single parameter to represent the total size of the
##' population -- over all compartments
##' @param infected name of a single parameter to represent the initial total size of
##' the infected population -- over all infected compartments
##' @family flexmodels
##' @export
initial_population = function(model, total, infected) {
    spec_check(
      introduced_version = '0.1.1',
      feature = 'Specification of total population and initial numbers of infected')
    model$initial_population = list(total = total, infected = infected)
    model
}

##' Add State Mappings
##'
##' Add regular expressions for identifying what states
##' are included in various stages of the eigenvector-based
##' method for constructing initial states
##' (\url{https://canmod.net/misc/flex_specs#v0.1.1})
##'
##' @param model \code{\link{flexmodel}} object
##' @param eigen_drop_pattern regular expression for identifying states
##' to be dropped before computing the eigenvector of the Jacobian
##' of the linearized model
##' @param infected_drop_pattern regular expression for identifying states
##' to be dropped from the eigenvector so that only 'infected' states
##' remain
##' @param initial_susceptible_pattern regular expression for
##' identifying states associated with susceptible classes
##'
##' @return \code{\link{flexmodel}} object
##'
##' @family flexmodels
##' @export
add_state_mappings = function(
    model,
    eigen_drop_pattern,
    infected_drop_pattern,
    initial_susceptible_pattern) {

    spec_check(
      introduced_version = '0.1.1',
      feature = 'Mapping between different vectors containing state information'
    )

    # TODO: pull outflow_to out of disease_free, because it
    # makes more sense as a general concept -- especially with
    # initial_susceptible_pattern -- and maybe in the future
    # it would be good to put things like vaccination categories
    # in there (OTOH vaccination categories are 'user-defined'
    # concepts whereas initial_susceptible_pattern and
    # infected_drop_pattern are 'general' concepts)
    #model$disease_free$state$patterns =
    model$initialization_mapping =
        list(eigen = eigen_drop_pattern,
             infected = infected_drop_pattern,
             susceptible = initial_susceptible_pattern)
    return(model)
}

#' Update Optimization Parameters
#'
#' \code{params ~ value}
#' \code{trans_param ~ value}
#' \code{trans_param ~ prior(hyperparameters ...)}
#' \code{trans_param ~ prior(hyperparameters ..., laplace = TRUE)}
#'
#' @export
update_opt_params = function(model, ...) {
  model$opt_params = lapply(list(...), parse_and_resolve_opt_form, model$params)
  model
}

#' @export
add_opt_params = function(model, ...) {
  # TODO: check for inconsistent specifications
  # (e.g. two different priors specified for beta0)
  model$opt_params = c(
    model$opt_params,
    lapply(list(...), parse_and_resolve_opt_form, model$params)
  )
  model
}

#' @export
initialize_opt_params = function(model) {
  # list(list( occurs because the result needs to
  # be a list of parsed formulas and each parsed
  # formula is itself a list
  list(list(
    param = list(
      param_nms = names(model$params),
      trans = rep('', length(model$params))
    ),
    prior = list(
      distr = 'flat',
      trans = '',
      reg_params = list(c(model$params))
    )
  ))
}

#' @export
initialize_opt_tv_params = function(model) {
  sc = model$timevar$piece_wise$schedule
  nms = unique(sc$Symbol)
  f = function(nm) {
    reg_params = (sc
      %>% filter(Symbol == nm)
      %>% getElement('init_tv_mult')
    )
    list(
      param = list(
        param_nms = nm,
        trans = ''
      ),
      prior = list(
        distr = 'flat',
        trans = '',
        reg_params = list(reg_params)
      )
    )
  }
  lapply(nms, f)
}

#' @export
update_opt_tv_params = function(
  model,
  tv_type = c('abs', 'rel_orig', 'rel_prev'),
  ...
) {
  tv_type = match.arg(tv_type, several.ok = TRUE)
  tvp = (model
      $  timevar
      $  piece_wise
      $  schedule
     %>% filter (Type %in% tv_type)
  )
  model$opt_tv_params = lapply(list(...), parse_and_resolve_opt_form, model$params)
  model
}

add_opt_tv_params = function(
  model,
  tv_type = c('abs', 'rel_orig', 'rel_prev'),
  ...
) {
  tv_type = match.arg(tv_type, several.ok = TRUE)
  tvp = (model
      $  timevar
      $  piece_wise
      $  schedule
     %>% filter (Type %in% tv_type)
  )
  model$opt_tv_params = c(
    model$opt_tv_params,
    lapply(list(...), parse_and_resolve_opt_form, model$params)
  )
  model
}

update_opt_vec = function(model, ...) {
  stop(
    "this is supposed to allow priors on (struc) vectors ",
    "and will wrap update_opt_params, but is not yet implemented. ",
    "this would allow multivariate priors, for example"
  )
}

# compute indices and pass them to the tmb/c++ side ---------------------

##' Update TMB Indices
##'
##' Add or update indices used to access appropriate values
##' during simulation and calibration using TMB
##'
##' \describe{
##'   \item{\code{make_ratemat_indices$from}}{}
##'   \item{\code{make_ratemat_indices$to}}{}
##'   \item{\code{make_ratemat_indices$count}}{}
##'   \item{\code{updateidx}}{indices into the \code{to}, \code{from}, and \code{count}
##'   vectors, identifying rates that need (or at least will) be updated at
##'   every simulation step}
##' }
##'
##' @param model compartmental model
##' @param another compartmental model with indices for TMB
##' @family flexmodels
##' @export
update_tmb_indices <- function(model) {

    # reduce rates so that there is only one rate
    # for each from-to pair
    model$rates = reduce_rates(model$rates)

    model$tmb_indices <- tmb_indices(model)
    return(model)
}

##' @inheritParams update_tmb_indices
##' @rdname update_tmb_indices
##' @export
add_tmb_indices = function(model) {
  stop("add_tmb_indices is no longer allowed. please use update_tmb_indices")
}

##' @family flexmodels
##' @export
tmb_indices <- function(model) {
    check_spec_ver_archived()
    if (length(model$rates) == 0L) {
      stop("no rates have been added to this model.\n",
           "please use one of add_rate, rep_rate, or vec_rate")
    }

    if (spec_ver_eq("0.1.0")) {
        # TODO: add model$factr here?
        sp <- c(model$state, model$params, model$sum_vector)
    } else {
        sp <- c(model$state, model$params)
    }

    indices = init_tmb_indices

    indices$make_ratemat_indices = ratemat_indices(model$rates, sp)

    if (spec_ver_gt("0.0.1")) {
        indices$par_accum_indices <-
            which(colnames(model$ratemat) %in% model$parallel_accumulators)
    }
    if (spec_ver_eq("0.0.2")) {
        indices$update_ratemat_indices <-
            ratemat_indices(state_dependent_rates(model), sp)
    }
    if (spec_ver_eq("0.0.3")) {
        indices$update_ratemat_indices <-
            ratemat_indices(time_varying_rates(model), sp)
    }
    if (spec_ver_gt("0.0.3")) {
        indices$updateidx <- which_time_varying_rates(model)
    }
    if (spec_ver_gt("0.0.6")) {
        indices$sum_indices = sum_indices(model$sums, model$state, model$params)
    }
    if (spec_ver_gt("0.1.0")) {
        if ((length(model$outflow) == 0L) & getOption("MP_warn_no_outflow")) {
          warning("model does not contain any outflow.\n",
                  "use add_outflow to balance inflows with outflows.\n",
                  "silence this warning with options(MP_warn_no_outflow = FALSE)")
        }

        indices$disease_free = disease_free_indices(model)
        indices$linearized_params = linearized_param_indices(model)
        indices$outflow = outflow_indices(model$outflow, model$ratemat)
        indices$linearized_outflow = outflow_indices(
            model$linearized_outflow, model$ratemat)
        indices$initialization_mapping = initialization_mapping_indices(model)
        indices$initial_population = initial_population_indices(model)
    }
    if (spec_ver_gt("0.1.1")) {
      indices$factr_indices = factr_indices(
        model$factrs,
        c(model$state, model$params, model$sum_vector)
      )
    }
    if (spec_ver_gt("0.1.2")) {
      indices$sim_report_expr_indices = sim_report_expr_indices(
        model$sim_report_exprs,
        initial_sim_report_names(model)
      )
      indices$lag_diff = lag_diff_indices(model)
      indices$conv = conv_indices(model)
      indices$observed = tmb_observed_data(model)
      indices$opt_params = tmb_opt_params(model)

      # indices into params vector and tv_mult vector
      # that identify parameters to be optimized
      opi = indices$opt_params$index_table$opt_param_id
      tvpi = indices$opt_params$index_tv_table$opt_tv_mult_id

      sc = model$timevar$piece_wise$schedule

      # MakeADFun map argument
      params_map = factor(
        rep(NA, length(model$params)),
        levels = seq_along(opi)
      )
      tv_mult_map = factor(
        rep(NA, nrow(model$timevar$piece_wise$schedule)),
        levels = seq_along(tvpi)
      )
      params_map[opi] = factor(levels(params_map))
      tv_mult_map[tvpi] = factor(levels(tv_mult_map))
      indices$ad_fun_map = list(
        params = params_map,
        tv_mult = tv_mult_map
      )
    }

    return(indices)
}

##' Make Objective Function with TMB
##'
##' Construct an objective function in TMB from a \code{flexmodel}
##' object. The behaviour of \code{tmb_fun} depends on \code{spec_version()}
##'
##' @param model object of class \code{flexmodel}
##' @family flexmodels
##' @importFrom TMB MakeADFun
##' @useDynLib McMasterPandemic
##' @export
tmb_fun <- function(model) {
    check_spec_ver_archived()
    DLL = getOption('MP_flex_spec_dll')

    if (getOption('MP_force_full_outflow')) {
      # reset outflow if it exists already
      if (length(model$outflow) > 0L) {
        model$outflow = list()
      }
      model = add_outflow(model)
    }
    if (getOption('MP_auto_tmb_index_update')) {
      model = update_tmb_indices(model)
    }

    unpack(model)
    unpack(tmb_indices)
    unpack(make_ratemat_indices)
    if (spec_ver_gt("0.0.3")) {
        unpack(timevar$piece_wise)
    }

    if (isTRUE(model$do_make_state)) {
      if (isTRUE(length(tmb_indices$disease_free$df_state_idx) == 0L)) {
        stop('cannot make the initial state because a disease-free state was not supplied. ',
             'either choose do_make_state = FALSE when initializing the model, ',
             'or use update_disease_free_state')
      }
    }


    ## make ad functions for different spec versions
    if (spec_ver_eq("0.0.1")) {
        dd <- MakeADFun(
            data = list(
                state = c(state),
                ratemat = ratemat,
                from = from,
                to = to,
                count = count,
                spi = spi,
                modifier = modifier
            ),
            parameters = list(params = c(unlist(params))),
            DLL = DLL
        )
    } else if (spec_ver_eq("0.0.2")) {
        up <- update_ratemat_indices

        ## TODO: spec problem -- numIters should have been kept in model
        numIters <- 3

        dd <- MakeADFun(
            data = list(
                state = c(state),
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
                numIterations = numIters
            ),
            parameters = list(params = c(unlist(params))),
            DLL = DLL
        )
    } else if (spec_ver_eq("0.0.4")) {
        dd <- MakeADFun(
            data = list(
                state = c(state),
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
                numIterations = iters
            ),
            parameters = list(params = c(unlist(params))),
            DLL = DLL
        )
    } else if (spec_ver_eq("0.0.5")) {
        dd <- MakeADFun(
            data = list(
                state = c(state),
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
                numIterations = iters
            ),
            parameters = list(params = c(unlist(params))),
            DLL = DLL
        )
    } else if (spec_ver_eq("0.0.6")) {
        dd <- MakeADFun(
            data = list(
                state = c(state),
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
                tv_mult = schedule$Value,
                tv_orig = schedule$Type == "rel_orig",
                par_accum_indices = par_accum_indices,
                do_hazard = do_hazard,
                numIterations = iters
            ),
            parameters = list(params = c(unlist(params))),
            DLL = DLL
        )
    } else if (spec_ver_eq("0.1.0")) {
        unpack(sum_indices)
        dd <- MakeADFun(
            data = list(
                state = c(state),
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
                tv_mult = schedule$Value,
                tv_orig = schedule$Type == "rel_orig",
                sumidx = sumidx,
                sumcount = unname(sumcount),
                summandidx = summandidx,
                par_accum_indices = par_accum_indices,
                do_hazard = do_hazard,
                numIterations = iters
            ),
            parameters = list(params = c(unlist(params))),
            DLL = DLL
        )
    } else if (spec_ver_eq("0.1.1")) {
        unpack(sum_indices)
        init_tv_mult = integer(0L)
        if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult
        dd <- MakeADFun(
            data = list(
                state = c(state),
                ratemat = ratemat,
                from = null_to_int0(from),
                to = null_to_int0(to),
                count = null_to_int0(count),
                spi = null_to_int0(spi),
                modifier = null_to_int0(modifier),
                updateidx = null_to_int0(c(updateidx)),
                breaks = null_to_int0(breaks),
                count_of_tv_at_breaks = null_to_int0(count_of_tv_at_breaks),
                tv_val = null_to_num0(schedule$tv_val),
                tv_spi = null_to_int0(schedule$tv_spi),
                tv_spi_unique = null_to_int0(sort(unique(schedule$tv_spi))),
                # tv_mult = schedule$Value,  # moved to parameter vector
                tv_orig = null_to_log0(schedule$Type == "rel_orig"),
                tv_abs = null_to_log0(schedule$Type == "abs"),

                sumidx = null_to_int0(sumidx),
                sumcount = null_to_int0(unname(sumcount)),
                summandidx = null_to_int0(summandidx),

                do_make_state = isTRUE(do_make_state),
                max_iters_eig_pow_meth = int0_to_0(null_to_0(max_iters_eig_pow_meth)),
                tol_eig_pow_meth = null_to_num0(tol_eig_pow_meth),

                outflow_row_count = null_to_int0(outflow$row_count),
                outflow_col_count = null_to_int0(outflow$col_count),
                outflow_rows = null_to_int0(outflow$rows),
                outflow_cols = null_to_int0(outflow$cols),

                linearized_outflow_row_count = null_to_int0(linearized_outflow$row_count),
                linearized_outflow_col_count = null_to_int0(linearized_outflow$col_count),
                linearized_outflow_rows = null_to_int0(linearized_outflow$rows),
                linearized_outflow_cols = null_to_int0(linearized_outflow$cols),

                lin_param_vals = null_to_num0(linearized_params$lin_param_vals),
                lin_param_count = null_to_int0(linearized_params$lin_param_count),
                lin_param_idx = null_to_int0(linearized_params$lin_param_idx),

                df_state_par_idx = null_to_int0(disease_free$df_state_par_idx),
                df_state_count = null_to_int0(disease_free$df_state_count),
                df_state_idx = null_to_int0(disease_free$df_state_idx),

                im_all_drop_eigen_idx = null_to_int0(initialization_mapping$all_drop_eigen_idx),
                im_eigen_drop_infected_idx = null_to_int0(initialization_mapping$eigen_drop_infected_idx),
                im_all_to_infected_idx = null_to_int0(initialization_mapping$all_to_infected_idx),
                im_susceptible_idx = null_to_int0(initialization_mapping$susceptible_idx),

                ip_total_idx = int0_to_0(null_to_0(initial_population$total_idx)),
                ip_infected_idx = int0_to_0(null_to_0(initial_population$infected_idx)),

                do_hazard = isTRUE(do_hazard),
                do_hazard_lin = isTRUE(do_hazard_lin),
                do_approx_hazard = isTRUE(do_approx_hazard),
                do_approx_hazard_lin = isTRUE(do_approx_hazard_lin),
                haz_eps = haz_eps,

                numIterations = int0_to_0(null_to_0(iters))
            ),
            parameters = list(params = c(unlist(params)),
                              tv_mult = init_tv_mult),
            DLL = DLL
        )
    } else if (spec_ver_eq("0.1.2")) {
      unpack(sum_indices)
      init_tv_mult = integer(0L)
      if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult
      dd <- MakeADFun(
        data = list(
          state = c(state),
          ratemat = ratemat,
          from = null_to_int0(from),
          to = null_to_int0(to),
          count = null_to_int0(count),
          spi = null_to_int0(spi),
          modifier = null_to_int0(modifier),
          updateidx = null_to_int0(c(updateidx)),
          breaks = null_to_int0(breaks),
          count_of_tv_at_breaks = null_to_int0(count_of_tv_at_breaks),
          tv_val = null_to_num0(schedule$tv_val),
          tv_spi = null_to_int0(schedule$tv_spi),
          tv_spi_unique = null_to_int0(sort(unique(schedule$tv_spi))),
          # tv_mult = schedule$Value,  # moved to parameter vector
          tv_orig = null_to_log0(schedule$Type == "rel_orig"),
          tv_abs = null_to_log0(schedule$Type == "abs"),

          sumidx = null_to_int0(sumidx),
          sumcount = null_to_int0(unname(sumcount)),
          summandidx = null_to_int0(summandidx),

          factr_spi = null_to_int0(factr_indices$spi_factr),
          factr_count = null_to_int0(factr_indices$count),
          factr_spi_compute = null_to_int0(factr_indices$spi),
          factr_modifier = null_to_int0(factr_indices$modifier),

          do_make_state = isTRUE(do_make_state),
          max_iters_eig_pow_meth = int0_to_0(null_to_0(max_iters_eig_pow_meth)),
          tol_eig_pow_meth = null_to_num0(tol_eig_pow_meth),

          outflow_row_count = null_to_int0(outflow$row_count),
          outflow_col_count = null_to_int0(outflow$col_count),
          outflow_rows = null_to_int0(outflow$rows),
          outflow_cols = null_to_int0(outflow$cols),

          linearized_outflow_row_count = null_to_int0(linearized_outflow$row_count),
          linearized_outflow_col_count = null_to_int0(linearized_outflow$col_count),
          linearized_outflow_rows = null_to_int0(linearized_outflow$rows),
          linearized_outflow_cols = null_to_int0(linearized_outflow$cols),

          lin_param_vals = null_to_num0(linearized_params$lin_param_vals),
          lin_param_count = null_to_int0(linearized_params$lin_param_count),
          lin_param_idx = null_to_int0(linearized_params$lin_param_idx),

          df_state_par_idx = null_to_int0(disease_free$df_state_par_idx),
          df_state_count = null_to_int0(disease_free$df_state_count),
          df_state_idx = null_to_int0(disease_free$df_state_idx),

          im_all_drop_eigen_idx = null_to_int0(initialization_mapping$all_drop_eigen_idx),
          im_eigen_drop_infected_idx = null_to_int0(initialization_mapping$eigen_drop_infected_idx),
          im_all_to_infected_idx = null_to_int0(initialization_mapping$all_to_infected_idx),
          im_susceptible_idx = null_to_int0(initialization_mapping$susceptible_idx),

          ip_total_idx = int0_to_0(null_to_0(initial_population$total_idx)),
          ip_infected_idx = int0_to_0(null_to_0(initial_population$infected_idx)),

          do_hazard = isTRUE(do_hazard),
          do_hazard_lin = isTRUE(do_hazard_lin),
          do_approx_hazard = isTRUE(do_approx_hazard),
          do_approx_hazard_lin = isTRUE(do_approx_hazard_lin),
          # haz_eps = haz_eps,

          sri_output = null_to_int0(sim_report_expr_indices$sri_output),
          sr_count = null_to_int0(sim_report_expr_indices$sr_count),
          sri = null_to_int0(sim_report_expr_indices$sri),
          sr_modifier = null_to_int0(sim_report_expr_indices$sr_modifier),

          lag_diff_sri = null_to_int0(lag_diff$sri),
          lag_diff_delay_n = null_to_int0(lag_diff$delay_n),

          conv_sri = null_to_int0(conv$sri),
          conv_c_prop_idx = null_to_int0(conv$c_prop_idx),
          conv_c_delay_cv_idx = null_to_int0(conv$c_delay_cv_idx),
          conv_c_delay_mean_idx = null_to_int0(conv$c_delay_mean_idx),
          conv_qmax = null_to_int0(conv$qmax),

          numIterations = int0_to_0(null_to_0(iters))
        ),
        parameters = list(params = c(unlist(params)),
                          tv_mult = init_tv_mult),
        DLL = DLL
      )

    } else if (spec_ver_gt("0.1.1")) {

      unpack(sum_indices)
      unpack(opt_params)

      # update parameter vectors -----------------
      init_tv_mult = integer(0L)
      if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult

      if(isTRUE(exists_opt_params(model))) {

        # tell tmb what parameters to put in the objective function
        map = ad_fun_map

        # pass transformed parameters to the objective function
        params = tmb_params_trans(model, vec_type = 'params')
        init_tv_mult = tmb_params_trans(model, vec_type = 'tv_mult')

      } else {
        map = list()
      }
      # -------------------------------------------

      if (!all(between(observed$time_step, 2, iters + 1))) {
        stop('observations outside of simulation range')
      }

      dd <- MakeADFun(
        data = list(
          state = c(state),
          ratemat = ratemat,
          from = null_to_int0(from),
          to = null_to_int0(to),
          count = null_to_int0(count),
          spi = null_to_int0(spi),
          modifier = null_to_int0(modifier),
          updateidx = null_to_int0(c(updateidx)),
          breaks = null_to_int0(breaks),
          count_of_tv_at_breaks = null_to_int0(count_of_tv_at_breaks),
          tv_val = null_to_num0(schedule$tv_val),
          tv_spi = null_to_int0(schedule$tv_spi),
          tv_spi_unique = null_to_int0(sort(unique(schedule$tv_spi))),
          # tv_mult = schedule$Value,  # moved to parameter vector
          tv_orig = null_to_log0(schedule$Type == "rel_orig"),
          tv_abs = null_to_log0(schedule$Type == "abs"),

          sumidx = null_to_int0(sumidx),
          sumcount = null_to_int0(unname(sumcount)),
          summandidx = null_to_int0(summandidx),

          factr_spi = null_to_int0(factr_indices$spi_factr),
          factr_count = null_to_int0(factr_indices$count),
          factr_spi_compute = null_to_int0(factr_indices$spi),
          factr_modifier = null_to_int0(factr_indices$modifier),

          do_make_state = isTRUE(do_make_state),
          max_iters_eig_pow_meth = int0_to_0(null_to_0(max_iters_eig_pow_meth)),
          tol_eig_pow_meth = null_to_num0(tol_eig_pow_meth),

          outflow_row_count = null_to_int0(outflow$row_count),
          outflow_col_count = null_to_int0(outflow$col_count),
          outflow_rows = null_to_int0(outflow$rows),
          outflow_cols = null_to_int0(outflow$cols),

          linearized_outflow_row_count = null_to_int0(linearized_outflow$row_count),
          linearized_outflow_col_count = null_to_int0(linearized_outflow$col_count),
          linearized_outflow_rows = null_to_int0(linearized_outflow$rows),
          linearized_outflow_cols = null_to_int0(linearized_outflow$cols),

          lin_param_vals = null_to_num0(linearized_params$lin_param_vals),
          lin_param_count = null_to_int0(linearized_params$lin_param_count),
          lin_param_idx = null_to_int0(linearized_params$lin_param_idx),

          df_state_par_idx = null_to_int0(disease_free$df_state_par_idx),
          df_state_count = null_to_int0(disease_free$df_state_count),
          df_state_idx = null_to_int0(disease_free$df_state_idx),

          im_all_drop_eigen_idx = null_to_int0(initialization_mapping$all_drop_eigen_idx),
          im_eigen_drop_infected_idx = null_to_int0(initialization_mapping$eigen_drop_infected_idx),
          im_all_to_infected_idx = null_to_int0(initialization_mapping$all_to_infected_idx),
          im_susceptible_idx = null_to_int0(initialization_mapping$susceptible_idx),

          ip_total_idx = int0_to_0(null_to_0(initial_population$total_idx)),
          ip_infected_idx = int0_to_0(null_to_0(initial_population$infected_idx)),

          do_hazard = isTRUE(do_hazard),
          do_hazard_lin = isTRUE(do_hazard_lin),
          do_approx_hazard = isTRUE(do_approx_hazard),
          do_approx_hazard_lin = isTRUE(do_approx_hazard_lin),
          # haz_eps = haz_eps,

          sri_output = null_to_int0(sim_report_expr_indices$sri_output),
          sr_count = null_to_int0(sim_report_expr_indices$sr_count),
          sri = null_to_int0(sim_report_expr_indices$sri),
          sr_modifier = null_to_int0(sim_report_expr_indices$sr_modifier),

          lag_diff_sri = null_to_int0(lag_diff$sri),
          lag_diff_delay_n = null_to_int0(lag_diff$delay_n),

          conv_sri = null_to_int0(conv$sri),
          conv_c_prop_idx = null_to_int0(conv$c_prop_idx),
          conv_c_delay_cv_idx = null_to_int0(conv$c_delay_cv_idx),
          conv_c_delay_mean_idx = null_to_int0(conv$c_delay_mean_idx),
          conv_qmax = null_to_int0(conv$qmax),

          obs_var_id = null_to_int0(observed$variable_id),
          obs_loss_id = null_to_int0(observed$loss_id),
          obs_loss_param_count = null_to_int0(observed$loss_param_count),
          obs_spi_loss_param = null_to_int0(observed$spi_loss_param),
          obs_time_step = null_to_int0(observed$time_step),
          obs_history_col_id = null_to_int0(observed$history_col_id),
          obs_value = observed$observed, # don't need to worry about missing values because they are omitted


          opt_param_id = null_to_int0(index_table$opt_param_id),
          opt_trans_id = null_to_int0(index_table$param_trans_id),
          opt_count_reg_params = null_to_int0(index_table$count_hyperparams),
          opt_reg_params = null_to_num0(hyperparameters),
          opt_reg_family_id = null_to_int0(index_table$prior_distr_id),

          opt_tv_param_id = null_to_int0(index_tv_table$opt_tv_mult_id),
          opt_tv_trans_id = null_to_int0(index_tv_table$param_trans_id),
          opt_tv_count_reg_params = null_to_int0(index_tv_table$count_hyperparams),
          opt_tv_reg_params = null_to_num0(hyperparameters_tv),
          opt_tv_reg_family_id = null_to_int0(index_tv_table$prior_distr_id),

          numIterations = int0_to_0(null_to_0(iters))
        ),
        parameters = list(params = c(unlist(params)),
                          tv_mult = null_to_num0(init_tv_mult)),
        map = map,
        DLL = DLL
      )
    } else {
        stop("Construction of TMB functions is not supported by your installation of MacPan")
    }
    return(dd)
}

##' @export
update_initial_state = function(model, silent = FALSE) {
  spec_check(
    introduced_version = '0.1.1',
    feature = 'Update initial state vector with C++ make_state'
  )
  if (!model$do_make_state) {
    if (!silent) warning('this model does not have the ability to update its own state')
    return(model)
  }
  model$state[] = initial_state_vector(model)
  model
}

# observed data ---------------------------------

#' @export
update_observed = function(model, data, error_dist) {
  stop("deprecated ... now in init_model and canned models")
  spec_check(
    introduced_version = '0.2.0',
    feature = 'comparison with observed data'
  )
  stopifnot(isTRUE(all.equal(c(names(data)), c("date", "var", "value"))))
  allvars = model$condensation_map[final_sim_report_names(model)]
  obsvars = unique(data$var)
  stopifnot(all(obsvars %in% allvars))
  model$observed$data = data
  model$observed$loss_params = data.frame(
    Parameter = "nb_disp", # only choice: dispersion
    Distribution = "nb",   # only choice: negative binomial
    Variable = obsvars
  )
  model$observed$timevar_error = (data.frame(
      Date = model$start_date
    )
    %>% cbind(model$observed$loss_params)
  )
  model
}

# time variation updates ------------------

#' @export
update_piece_wise = function(model, params_timevar) {
  spec_check(
    introduced_version = '0.0.4',
    feature = 'piece-wise time variation of parameters'
  )

  schedule <- (params_timevar
               %>% mutate(Date = as.Date(Date))
               %>% mutate(breaks = (Date
                  %>% difftime(model$start_date, units = 'days')
                  %>% as.integer
               ))
               %>% mutate(tv_spi = find_vec_indices(Symbol, c(model$state, model$params)))
               %>% arrange(breaks, tv_spi)
  )

    if(spec_ver_gt("0.1.0")) {
      schedule = (schedule
                  %>% mutate(init_tv_mult = replace(Value,
                                                    which(is.na(Value)),
                                                    1))
                  %>% mutate(last_tv_mult = init_tv_mult)
      )
    }

    count_of_tv_at_breaks <- c(table(schedule$breaks))

    schedule$tv_val <- NA
    new_param <- TRUE
    ns <- nrow(schedule)
    for (i in 1:ns) {
      if (new_param | schedule$Type[i] == "rel_orig") {
        old_val <- model$params[schedule$Symbol[i]]
      } else if (schedule$Type[i] == "rel_prev") {
        old_val <- schedule$tv_val[i - 1]
      } else if (schedule$Type[i] == "abs") {
        old_val <- 1
      } else {
        stop(
          "Unrecognized break-point type.\n",
          "Only rel_orig, rel_prev, and abs are allowed"
        )
      }
      schedule$tv_val[i] <- old_val * schedule$Value[i]
      new_param <- schedule$Symbol[i] != schedule$Symbol[min(ns, i + i)]
    }

    attr(model$params, "tv_param_indices") <- setNames(find_vec_indices(
      unique(schedule$Symbol),
      model$params
    ), unique(schedule$Symbol))


    model$timevar$piece_wise <- list(
      ## schedule includes tv_spi and tv_val as in spec
      schedule = schedule,
      breaks = as.integer(names(count_of_tv_at_breaks)),
      count_of_tv_at_breaks = unname(count_of_tv_at_breaks)
    )

    if (getOption("MP_warn_bad_breaks")) {
      good_dates = between(
        model$timevar$piece_wise$schedule$Date,
        model$start_date,
        model$end_date - 1
      )
      if (!all(good_dates)) {
        warning("some time-varying parameters will not change at every\n",
                "specified date because some parameter changes only take\n",
                "effect outside of the simulation dates.\n",
                "to silence this warning use:\n",
                "options(MP_warn_bad_breaks = FALSE)")
      }
    }
    model
}

#' @export
initialize_piece_wise = function(model) {
  model$timevar = list(
    piece_wise = list(
      breaks = integer(0L),
      count_of_tv_at_breaks = integer(0L),
      schedule = data.frame(
        Date = as.Date(numeric(0L), origin = "1970-01-01"),
        Symbol = character(0L),
        Type = character(0L),
        Value = numeric(0L),
        tv_spi = integer(0L),
        tv_val = numeric(0L)
      )
    )
  )
  model
}

# regenerate model --------------------------------

if(FALSE) {

model = yukon_model

model %>% names

arg_names = function(f, include_dots = FALSE, warn_dots = TRUE) {
  nms = names(formals(f))
  if (include_dots) return(nms)
  is_dots = nms == "..."
  if (warn_dots & any(is_dots)) warning("removing dots")
  nms[!is_dots]
}

arg_names_classify = function(f, nms) {
  arg_nms = arg_names(f, TRUE)
  args_in_nms = arg_nms %in% nms
  nms_in_args = nms %in% arg_nms
  list(
    args_in_nms = arg_nms[args_in_nms],
    args_not_in_nms = arg_nms[!args_in_nms],
    nms_in_args = nms[nms_in_args],
    nms_not_in_args = nms[!nms_in_args]
  )
}

names(model)
params_timevar = model$timevar$piece_wise$schedule[c("Date", "Symbol", "Value", "Type")]
data = model$observed$data
arg_names_classify(init_model, names(model))
arg_names_classify(add_rate, names(model$rates[[1]]))
arg_names_classify(add_state_param_sum, names(model$sums[[1]]))
arg_names_classify(add_factr, names(model))
arg_names_classify(add_sim_report_expr, names(model))
arg_names_classify(add_lag_diff, names(model))
arg_names_classify(add_conv, names(model))
arg_names_classify(add_linearized_outflow, names(model))
arg_names_classify(add_outflow, names(model))
arg_names_classify(update_linearized_params, names(model))
arg_names_classify(update_disease_free_state, names(model))
arg_names_classify(initial_population, names(model))
arg_names_classify(add_state_mappings, names(model))
arg_names_classify(update_opt_params, names(model))
arg_names_classify(update_opt_tv_params, names(model))
}

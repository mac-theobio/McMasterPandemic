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
##' @param do_approx_hazard_lin like \code{do_approx_hazard} but for
##' the linearized model that is used to construct the initial
##' state (experimental)
##' @param do_make_state should state be remade on the c++ size?
##' (https://canmod.net/misc/flex_specs#v0.1.1) -- only used
##' if \code{spec_ver_gt('0.1.0')}
##' @param do_sim_constraint should simulated values be smoothly
##' constrained to be above \code{sim_lower_bound} when computing
##' negative binomial (maybe others in the future?) loss functions?
##' the smooth constraint function is
##' \eqn{y = x + \epsilon * \exp(-x / \epsilon)}, where \eqn{\epsilon}
##' is \code{sim_lower_bound} and \eqn{x} is the simulated value.
##' @param sim_lower_bound optional lower bound on the simulated values
##' when computing negative binomial loss functions (only applicable
##' when \code{do_sim_constraint} is \code{TRUE}.
##' @param max_iters_eig_pow_meth maximum number of iterations
##' to use in computing the eigenvector for initial state
##' construction
##' @param tol_eig_pow_meth tolerance for determining convergence
##' of the power method used in initial state construction
##' @param data optional observed data frame in long format to
##' compare with simulated trajectories. must have the following
##' columns: \code{date}, \code{var}, \code{value}. (currently this is not working)
##' @family flexmodels
##' @return flexmodel object representing a compartmental model
##' @importFrom lubridate Date
##' @export
flexmodel <- function(params, state = NULL,
       start_date = NULL, end_date = NULL,
       params_timevar = NULL,
       do_hazard = getOption("MP_default_do_hazard"),
       do_make_state = getOption("MP_default_do_make_state"),
       do_hazard_lin = getOption("MP_default_do_hazard_lin"),
       do_approx_hazard = getOption("MP_default_do_approx_hazard"),
       do_approx_hazard_lin = getOption("MP_default_do_approx_hazard_lin"),
       do_sim_constraint = getOption("MP_default_do_sim_constraint"),
       sim_lower_bound = getOption("MP_default_sim_lower_bound"),
       max_iters_eig_pow_meth = 8000,
       tol_eig_pow_meth = 1e-6,
       data = NULL, ...) {
    if(!is.null(data)) {
      stop(
        "currently the data argument is not working.\n",
        "please use update_observed instead"
      )
    }
    check_spec_ver_archived()
    name_regex = wrap_exact(getOption("MP_name_search_regex"))
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
      op = options(MP_rexp_steps_default = 1)
      state = make_state(params = params)
      options(op)
      state[] = 0
    } else if (is.character(state)) {
      state = const_named_vector(state, 0)
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

    params = unlist_params(params)

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
        model$iters <- compute_num_iters(model)
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
          model = update_piece_wise(model, params_timevar, regenerate_rates = FALSE)
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
        model$haz_eps = 1e-6
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
    model$pow = init_pow
    model$pow_vector = init_pow_vector

    # condensation -- spec_ver_gt("0.1.2)
    model$condensation = init_condensation


    if (spec_ver_gt("0.1.2") & !is.null(data)) {
      model = update_observed(model, data)
    } else {
      model$observed = init_observed
    }

    if (spec_ver_gt('0.1.2')) {
      model$do_sim_constraint <- do_sim_constraint
      model$sim_lower_bound <- sim_lower_bound
    }

    model$opt_params = list()
    model$opt_tv_params = list(
      abs = list(),
      rel_orig = list(),
      rel_prev = list()
    )
    model$condensation_map = init_condensation_map
    model$no_condensation = TRUE

    if (FALSE & spec_ver_gt('0.1.2')) {
      model$opt_params = initialize_opt_params(model)
      model$opt_tv_params = initialize_opt_tv_params(model)
    }

    model$tmb_indices <- init_tmb_indices
    class(model) = 'flexmodel'

    model
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

## Define Rate for Single Element of Rate Matrix
##
## @param from from state
## @param to to state
## @param formula one-sided formula defining the rate with reference
## to the parameters and state variables
## @param state state_pansim object
## @param params param_pansim object
## @param sums vector of sums of state variables and parameters
## @param factrs vector of factrs ...
## @param pows vector of pows ...
## @param ratemat rate matrix
## @importFrom dplyr bind_rows
## @family flexmodel_definition_functions

rate <- function(from, to, formula, state, params, sums, factrs, pows, ratemat) {
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
            c(state, params, sums, factrs, pows))

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

##' Add Rate
##'
##' Define how the rate of flow from one compartment to another
##' depends on parameters, state variables, .
##'
##' @param model \code{\link{flexmodel}} object
##' @param from Name of state from which flow is coming
##' @param to Name of state to which flow is going
##' @param formula Model formula defining dependence of the rate on
##' existing variables. See \code{\link{avail_for_rate}} for a function
##' that will print out the names of variables that are available for
##' use in these formulas.
##' @return updated \code{\link{flexmodel}} with an additional
##' non-zero rate matrix
##' element specified
##' @family flexmodel_definition_functions
##' @family rate_functions
##' @export
add_rate <- function(model, from, to, formula) {
    unpack(model)
    added_rate <- (from
      %>% rate(to, formula, state, params, sum_vector, factr_vector, pow_vector, ratemat)
      %>% list
      %>% setNames(paste(from, to, sep = "_to_"))
    )
    model$rates <- c(model$rates, added_rate)
    return(model)
}

#' Repeat a Rate for Several Rate Matrix Elements
#'
#' @param model \code{\link{flexmodel}} object
#' @param from character vector defining states from which flow is coming
#' @param to character vector defining states from which flow is going
#' @param formula Model formula defining dependence of the repeated rate on
#' existing variables. See \code{\link{avail_for_rate}} for a function
#' that will print out the names of variables that are available for
#' use in these formulas.
#' @param mapping experimental -- please choose default for now
#' @family flexmodel_definition_functions
#' @family rate_functions
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
        MoreArgs = nlist(formula, state, params, sums,
                         factrs = factr_vector, pows = pow_vector,
                         ratemat),
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    nms = mapply(paste, from, to, MoreArgs = list(sep = "_to_"))
    model$rates <- c(rates, setNames(lst, nms))

    return(model)
}

#' Specify Vector of Rates
#'
#' @param model \code{\link{flexmodel}} object
#' @param from character vector defining states from which flow is coming
#' @param to character vector defining states from which flow is going
#' @param formula \code{\link{struc-class}} object defining a vector of flows
#' for each \code{from-to} pair. See \code{\link{avail_for_rate}} for a
#' function that will print out the names of variables that are available for
#' use in these \code{\link{struc-class}} objects.
#' @param mapping experimental -- please choose default for now
#' @family flexmodel_definition_functions
#' @family rate_functions
#' @export
vec_rate = function(model, from, to, formula) {

    unpack(model)
    #check_from_to(from, to, names(state))

    indices = pwise(from, to, ratemat)

    from = rownames(ratemat)[indices[,'from_pos']]
    to = colnames(ratemat)[indices[,'to_pos']]

    lst = mapply(rate, from, to, as.character(formula),
                 MoreArgs = nlist(state, params, sums, factrs = factr_vector, pows = pow_vector, ratemat),
                 SIMPLIFY = FALSE, USE.NAMES = FALSE)
    nms = mapply(paste, from, to, MoreArgs = list(sep = "_to_"))
    model$rates <- c(rates, setNames(lst, nms))

    return(model)
}

#' Specify Matrix of Rates
#'
#' Not implemented
#'
#' @family flexmodel_definition_functions
#' @family rate_functions
#' @export
mat_rate = function() {
    stop("\nrate specification with matrices is ",
         "coming sometime in the future ... maybe\n",
         "in the meantime you can specify vector-valued rates with vec_rate\n",
         "see this document for potentially more information on priorities:\n",
         options("MP_flex_spec_doc_site")[[1]])
}



# factr and associated functions ----------------------

factr <- function(factr_nm, formula, state, params, sums, factrs, ratemat) {
  ## TODO: test for formula structure
  stopifnot(
    (  inherits(formula, "formula")
       | is.character(formula)
       | is(formula, 'struc')
    )
  )

  product_list <- function(x) {
    #factor_table = McMasterPandemic::factor_table
    #find_vec_indices = McMasterPandemic::find_vec_indices
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
#' or 1-by-1 \code{\link{struc-class}} object describing the
#' intermediate factor. See \code{\link{avail_for_factr}} for
#' a function that will return the names of all variables that are available
#' for use in these formulas and \code{\link{struc-class}} objects.
#'
#' @seealso See \code{\link{vec_factr}} to add more than
#' one intermediate factor at the same time.
#'
#' @return updated \code{\link{flexmodel}} object
#' @family flexmodel_definition_functions
#' @export
add_factr <- function(model, factr_nm, formula) {
  state <- params <- sum_vector <- factr_vector <- state <- params <-
    sum_vector <- factr_vector <- pow_vector <- ratemat <- NULL
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
#' @param formula \code{\link{struc-class}} object describing the
#' vector of intermediate factors. See \code{\link{avail_for_factr}} for
#' a function that will return the names of all variables that are available
#' for use in these \code{\link{struc-class}} objects.
#'
#' @seealso See \code{\link{add_factr}} to add a single
#' scalar-valued intermediate factor and \code{\link{vec}}
#' to create a vector-valued \code{\link{struc}} object.
#'
#' @return updated \code{\link{flexmodel}} object
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
      stop("sums cannot have the same name as state variables ",
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
#' @param model \code{\link{flexmodel}} object
#' @param sum_name name of sum of state variables and parameters
#' @param summands character vector of regular expressions for identifying
#' state variables and parameters to sum together. See \code{\link{avail_for_sum}}
#' for a function that will return all names that these regular expressions
#' will search through.
#'
#' @return updated \code{\link{flexmodel}}
#'
#' @family flexmodel_definition_functions
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

pow = function(pow_nms, pow_arg1_nms, pow_arg2_nms, pow_const_nms,
               state, params, sum_vector, factr_vector, pow_vector) {
  avail_vec = c(state, params, sum_vector, factr_vector, pow_vector)
  arg1_idx = find_vec_indices(pow_arg1_nms, avail_vec)
  arg2_idx = find_vec_indices(pow_arg2_nms, avail_vec)
  const_idx = find_vec_indices(pow_const_nms, avail_vec)
  initial_value = setNames(
    avail_vec[const_idx] * (avail_vec[arg1_idx]^avail_vec[arg2_idx]),
    pow_nms
  )
  data.frame(pow_nms, pow_arg1_nms, pow_arg2_nms, pow_const_nms, initial_value)
}

#' Add Power Law
#'
#' Compute and save the intermediate result of applying a
#' power law to the state variables, parameters, sums of these
#' things, and factr expressions of them as well.
#'
#' @param model \code{\link{flexmodel}} object
#' @param pow_nms name of the resulting power law
#' @param pow_arg1_nms names of the variables used as the bases
#' @param pow_arg2_nms names of the variables used as the exponents
#' @param pow_const_nms names of the variables used as the constant
#'
#' @export
add_pow = function(model, pow_nms, pow_arg1_nms, pow_arg2_nms, pow_const_nms) {
  model$pow = rbind(model$pow, pow(
    pow_nms, pow_arg1_nms, pow_arg2_nms, pow_const_nms,
    model$state, model$params, model$sum_vector,
    model$factr_vector, model$pow_vector
  ))
  model$pow_vector = setNames(model$pow$initial_value, model$pow$pow_nms)
  model
}

# condensation -----------------------------------

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

#' Add Expression to the Simulation History
#'
#' Create a new variable in the \code{\link{simulation_history}}
#' by taking sums and products of existing variables (and their
#' complements and inverses). These new variables can be compared
#' with observed data streams.
#'
#' @param model \code{\link{flexmodel}} object
#' @param expr_nm Name of the new variable
#' @param formula formula describing the expression. See
#' \code{\link{avail_for_expr}} for
#' a function that will return the names of all variables that are available
#' for use in these formulas.
#'
#' @return updated \code{\link{flexmodel}} object
#'
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

#' Add Variable by Lag Differencing
#'
#' Add a variable that will be available for comparisons with observed data
#' that is created by processing an existing variable in the simulation
#' history by lagged differencing.
#'
#' @param model \code{\link{flexmodel}} object
#' @param var_pattern regular expression used to identify variables for
#' differencing. See \code{\link{avail_for_lag}} for a function that will
#' return the names of all variables that are available for differencing.
#' @param delay_n Delay in days for determining the lag in the differences.
#' @param lag_dates Dates between which differences should be taken.
#' @param input_names names of variables to be differenced
#' @param output_names names of the result of differencing
#'
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
  if (any(duplicated(unlist(lapply(model$lag_diff, getElement, "output_names"))))) {
    stop('the same variable is being differenced with the same lag, which is not currently allowed')
  }
  update_tmb_indices(model)
}

#' @rdname add_lag_diff
#' @export
add_lag_diff_uneven = function(model, input_names, output_names, lag_dates) {
  spec_check(
    introduced_version = '0.2.1',
    feature = 'uneven lagged differencing',
    exception_type = 'warning'
  )
  if (spec_ver_gt('0.2.0')) {
    stopifnot(input_names %in% intermediate_sim_report_names(model))
    stopifnot(all(lag_dates %in% simulation_dates(model)))
    model$lag_diff_uneven = c(
      model$lag_diff_uneven,
      list(nlist(input_names, output_names, lag_dates))
    )
  }
  model
}

#' Add Variable by Convolution
#'
#' Add a variable that will be available for comparisons with observed data
#' that is created by processing an existing variable in the simulation
#' history by convolving it with a gamma density.
#'
#' @param model \code{\link{flexmodel}} object
#' @param var_pattern regular expression used to identify variables for
#' convolutions. See \code{\link{avail_for_conv}} for a function that will
#' return the names of all variables that are available for convolution.
#' @param c_prop name of the parameter in the \code{params} vector associated
#' with the proportion of individuals in the original state variable that are
#' represented in the convolved variable
#' @param c_delay_cv name of the parameter in the \code{params} vector
#' associated with the coefficient of variation of the gamma density
#' @param c_delay_mean name of the parameter in the \code{params} vector
#' associated with the mean of the gamma density
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
  if (any(duplicated(unlist(lapply(model$conv, getElement, "output_names"))))) {
    stop('the same variable is being convolved twice, which is not currently allowed')
  }
  update_tmb_indices(model)
}

#' Specify and Name Variables for Condensation
#'
#' The condensed simulation is a (possibly renamed) subset
#' of the variables in the simulation history. This function
#' allows one to define this subset and how it is renamed,
#' by creating a map with the following structure:
#' \code{c(orig_var_i = "cond_var_1", ..., orig_var_j = "cond_var_n")}.
#'
#' @param model \code{\link{flexmodel}} object
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
  model$no_condensation = FALSE
  model
}


# parallel accumulators (deprecated -- use outflow instead) ---------------------

##' Add Parallel Accumulators
##'
##' [deprecated] Add parallel accumulators to a compartmental model.
##'
##' @param model \code{\link{flexmodel}} object
##' @param state_patterns regular expressions for identifying states as
##' parallel accumulators
##' @return updated \code{\link{flexmodel}} with parallel
##' accumulators specified
##' @family flexmodel_definition_functions
##' @export
add_parallel_accumulators <- function(model, state_patterns) {
    model$parallel_accumulators <- parallel_accumulators(model, state_patterns)
    return(model)
}

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
##' @return updated \code{\link{flexmodel}} object
##'
##' @family flexmodel_definition_functions
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

##' Update Linearized Model Parameters
##'
##' Specify how to update model parameters for use with
##' linearized models during initial state vector construction.
##'
##' @param model \code{\link{flexmodel}} object
##' @param param_pattern regular expression identifying parameters
##' that require updating before they can be used in linearized
##' model simulations
##' @param value numeric value required to update
##'
##' @return updated \code{\link{flexmodel}} object
##'
##' @family flexmodel_definition_functions
##' @export
update_linearized_params = function(model, param_pattern, value) {
    model$linearized_params = c(
        model$linearized_params,
        list(linearized_params(model, param_pattern, value)))
    return(model)
}

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

##' Update Disease-Free State
##'
##' @param model \code{\link{flexmodel}} object
##' @param state_pattern regular expression for identifying state variables
##' to be updated when constructing a disease-free state
##' @param param_pattern regular expression for identifying parameters
##' to use as disease-free state variables
##'
##' @return updated \code{\link{flexmodel}} object
##'
##' @family flexmodel_definition_functions
##' @export
update_disease_free_state = function(model, state_pattern, param_pattern) {
    model$disease_free = c(
        model$disease_free,
        list(disease_free_state(model, state_pattern, param_pattern)))
    return(model)
}

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
##' @param model \code{\link{flexmodel}} object
##' @param total name of a single parameter to represent the total size of the
##' population -- over all compartments
##' @param infected name of a single parameter to represent the initial total size of
##' the infected population -- over all infected compartments
##'
##' @return updated \code{\link{flexmodel}} object
##'
##' @family flexmodel_definition_functions
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
##' @return updated \code{\link{flexmodel}} object
##'
##' @family flexmodel_definition_functions
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



#' Optimization Parameters
#'
#' Add or update parameters to be optimized/calibrated in a
#' \code{flexmodel} object.
#'
#' @section Formula Syntax:
#' The left-hand-side of the formula describes the parameters and how they
#' should be transformed, and the right-hand-side describes the associated
#' prior distribution.
#'
#' The simplest approach is to specify one parameter at a time with the
#' following syntax:
#'
#' \code{trans_param ~ trans_prior(hyperparameters ...)}
#'
#' On the left-hand-side the optional transformation, \code{trans}, is
#' separated by an underscore from the parameter name, \code{param}. For
#' example one might optimize the transmission rate, \code{beta}, on the
#' log scale by specifying the left-hand-side as \code{log_beta}. To
#' optimize \code{beta} on the untransformed scale the left-hand-side
#' would be simply be \code{beta}. The currently available transformations
#' are \code{log}, \code{log10}, \code{logit}, \code{cloglog}, \code{inverse}.
#'
#' On the right-hand-side the optional transformation, \code{trans}, is
#' separated from the name of the prior family by an underscore. This
#' transformation defines the scale on which the prior distribution is over.
#' Currently the transformation scale passed to the objective function
#' (on the left-hand-side) must match the transformation scale of the
#' priot (on the right-hand-side) -- perhaps this restriction will be
#' lifted one day. The currently available prior families are \code{flat}
#' (for no regularization) and \code{normal}. These families are
#' expressed as functions of hyperparameters where the first argument
#' is a location parameter, which also defines the starting point of the
#' optimizer.
#'
#' \describe{
#'   \item{\code{flat}}{
#'     \describe{
#'       \item{initial}{the 'peak' of a flat function -- important because it sets the initial value of the optimizer.}
#'     }
#'   }
#'   \item{\code{normal}}{
#'       \describe{
#'         \item{mean}{mean of the normal distribution, also used as the initial value of the optimizer}
#'         \item{standard deviation}{standard deviation of the normal distribution}
#'       }
#'     }
#' }
#'
#' @param model \code{\link{flexmodel}} object
#' @param tv_type type of time-variation for time-varying parameters,
#' which can be one of \code{'abs', 'rel_orig', 'rel_prev',
#' 'rel_orig_logit', 'rel_orig_prev'}
#' @param ... a list of formulas for describing what parameters should
#' be optimized, whether/how they should be transformed before being
#' passed to the objective function, and what prior distribution (or
#' regularization function should be used -- see the section on
#' the formula syntax for more details.
#'
#' @return updated \code{\link{flexmodel}} object
#'
#' @export
add_opt_params = function(model, ...) {
  if (!inherits(model, "flexmodel_to_calibrate")) {
    stop("\nplease add observed data for model fitting, ",
         "\nbefore specifying parameters to be optimized"
    )
  }
  l = force(list(...))
  model$opt_params = c(
    model$opt_params,
    lapply(list(...), parse_and_resolve_opt_form, model$params)
  )
  update_tmb_indices(model)
}

#' @rdname add_opt_params
#' @export
update_opt_params = function(model, ...) {
  if (!inherits(model, "flexmodel_to_calibrate")) {
    stop("\nplease add observed data for model fitting, ",
         "\nbefore specifying parameters to be optimized"
    )
  }
  model$opt_params = lapply(list(...), parse_and_resolve_opt_form, model$params)
  update_tmb_indices(model)
}

#' @rdname add_opt_params
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

#' @rdname add_opt_params
#' @export
initialize_opt_tv_params = function(model) {
  # FIXME: doesn't seem necessary
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

#' @rdname add_opt_params
#' @export
update_opt_tv_params = function(
  model,
  tv_type = valid_tv_types,
  ...
) {
  if (!inherits(model, "flexmodel_to_calibrate")) {
    stop("\nplease add observed data for model fitting, ",
         "\nbefore specifying parameters to be optimized"
    )
  }
  tv_type = match.arg(tv_type, several.ok = FALSE)
  # tvp = (model
  #     $  timevar
  #     $  piece_wise
  #     $  schedule
  #    %>% filter (Type %in% tv_type)
  # )
  model$opt_tv_params[[tv_type]] = lapply(
    list(...),
    parse_and_resolve_opt_form, model$params
  )
  model
}

#' @rdname add_opt_params
#' @export
add_opt_tv_params = function(
  model,
  tv_type = valid_tv_types,
  ...
) {
  if (!inherits(model, "flexmodel_to_calibrate")) {
    stop("\nplease add observed data for model fitting, ",
         "\nbefore specifying parameters to be optimized"
    )
  }
  tv_type = match.arg(tv_type, several.ok = FALSE)
  # tvp = (model
  #     $  timevar
  #     $  piece_wise
  #     $  schedule
  #    %>% filter (Type %in% tv_type)
  # )
  model$opt_tv_params[[tv_type]] = c(
    model$opt_tv_params[[tv_type]],
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

# refine_initial_values = function(model)

##' Update Random Effect Parameters
##'
##' Specify parameters as random effects
##'
##' @param model \code{flexmodel_to_calibrate} object
##' @param ... formulas for specifying parameters as random effects -- see
##' \code{\link{update_opt_params}} for more detail
##'
##' @export
update_ranef_params = function(model, ...) {
  if (!inherits(model, "flexmodel_to_calibrate")) {
    stop("\nplease add observed data for model fitting, ",
         "\nbefore specifying parameters as random effects"
    )
  }
  model$ranef_params = lapply(list(...), parse_and_resolve_opt_form, model$params)
  model
}

#' Extend End Date
#'
#' Extend the final simulation date a number of days
#' in the future.
#'
#' @param model \code{\link{flexmodel}} object
#' @param days_to_extend number of days to extend the end date
#' @importFrom lubridate days
#' @export
extend_end_date = function(model, days_to_extend) {
  UseMethod('extend_end_date')
}

#' @exportS3Method
extend_end_date.flexmodel = function(model, days_to_extend) {
  model$end_date = model$end_date + days(days_to_extend)
  model$iters = compute_num_iters(model)
  update_tmb_indices(model)
}

#' @exportS3Method
extend_end_date.flexmodel_calibrated = function(model, days_to_extend) {
  model$model_to_calibrate = extend_end_date(
    model$model_to_calibrate,
    days_to_extend
  )
  NextMethod("extend_end_date")
}

#' Update Simulation Bounds
#'
#' @param model \code{\link{flexmodel}} object
#' @param start_date optional new start date for simulations
#' @param end_date optional new end date for simulations
#'
#' @export
update_simulation_bounds = function(model, start_date = NULL, end_date = NULL) {
  if (inherits(model, 'flexmodel_calibrated')) {
    stop("it is not currently allowed to update the simulation bounds on a calibrated model. please see ?extend_end_date for a possible solution")
  }
  if (!is.null(start_date)) {
    model$start_date = as.Date(start_date)
  }
  if (!is.null(end_date)) {
    model$end_date = as.Date(end_date)
  }
  model$iters = compute_num_iters(model)
  # TODO: check time-variation breakpoints? maybe updating the indices will be enough
  update_tmb_indices(model)
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
##' @param model \code{\link{flexmodel}} object
##' @export
update_tmb_indices = function(model) {
  if ((!spec_ver_eq(model$spec_ver)) & getOption("MP_no_indices_w_mismatched_specs")) {
    stop(
      "\nthis model was initialized under spec version ",
      model$spec_ver,
      ",\nbut spec version ",
      getOption("MP_flex_spec_version"),
      " is currently running\n(see ", getOption("MP_flex_spec_doc_site"),
      " for context).",
      "\nyou can set options(MP_no_indices_w_mismatched_specs = FALSE),",
      "\nbut this can be dangerous and lead to confusing errors.",
      "\nit is better to create the model again in the current environment.",
      "\nif you know what you are doing you can change the spec with",
      "\noptions(MP_flex_spec_version = ",
      sQuote(model$spec_ver), ") or set_spec_version(",
      sQuote(model$spec_ver), ", ...)"
    )
  }
  UseMethod("update_tmb_indices")
}

##' @exportS3Method
update_tmb_indices.flexmodel_to_calibrate = function(model) {
  data_vars = unique(model$observed$data$var)
  loss_vars = unique(model$observed$loss_params$Variable)
  if (!all(data_vars %in% loss_vars)) {
    # protects against segfaults!
    stop(
      "\nall variables in the observed data must be associated",
      "\nwith observation error. please use update_error_dist and",
      "\nadd_error_dist."
    )
  }
  NextMethod('update_tmb_indices')
}

##' @exportS3Method
update_tmb_indices.flexmodel <- function(model) {

    if (model$no_condensation) {
      model = update_full_condensation_map(model)
    }

    # reduce rates so that there is only one rate
    # for each from-to pair
    model$rates = reduce_rates(model$rates)

    model$tmb_indices <- tmb_indices(model)
    return(model)
}

##' @rdname update_tmb_indices
##' @export
add_tmb_indices = function(model) {
  stop("add_tmb_indices is no longer allowed. please use update_tmb_indices")
}

##' @family flexmodel_definition_functions
##' @rdname update_tmb_indices
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
        if ((length(model$outflow) == 0L) & getOption("MP_auto_outflow")) {
          model = add_outflow(model)
        }
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
      indices$pow_indices = pow_indices(
        model$pow,
        c(model$state, model$params, model$sum_vector, model$factr_vector, model$pow_vector)
      )
    }
    if (spec_ver_gt("0.1.2")) {
      indices$sim_report_expr_indices = sim_report_expr_indices(
        model$sim_report_exprs,
        initial_sim_report_names(model)
      )
      if (spec_ver_gt("0.2.0")) {
        indices$lag_diff = lag_diff_uneven_indices(model)
      } else {
        indices$lag_diff = lag_diff_indices(model)
      }
      indices$conv = conv_indices(model)
      indices$observed = tmb_observed_data(model)
      indices$opt_params = tmb_opt_params(model)

      # indices into params vector and tv_mult vector
      # that identify parameters to be optimized
      opi = indices$opt_params$index_table$opt_param_id
      tvpi = indices$opt_params$index_tv_table$opt_tv_mult_id

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
##' @param model \code{\link{flexmodel}} object
##' @importFrom TMB MakeADFun
##' @useDynLib McMasterPandemic
##' @export
tmb_fun = function(model) {
  do.call(MakeADFun, tmb_fun_args(model))
}

##' @rdname tmb_fun
##' @export
tmb_fun_args = function(model) {
  spec_with_underscores = gsub(
    "\\.", "_",
    getOption("MP_flex_spec_version")
  )
  class(model) = c(
    class(model),
    'spec_ver' %_% spec_with_underscores
  )
  tmb_fun_args_by_spec(model)
}

tmb_fun_args_by_spec = function(model) {
  check_spec_ver_archived()
  DLL = getOption('MP_flex_spec_dll')
  model = update_model_before_tmb_args(model)

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
  UseMethod("tmb_fun_args_by_spec")
}

update_model_before_tmb_args = function(model) {
    if (getOption('MP_auto_outflow') & spec_ver_gt("0.1.0")) {
      if (length(model$outflow) == 0L) {
        model = add_outflow(model)
      }
    }
    if (model$no_condensation) {
      model = update_full_condensation_map(model)
    }
    if (getOption('MP_auto_tmb_index_update')) {
      model = update_tmb_indices(model)
    }
    model
}

tmb_fun_args_by_spec.spec_ver_0_0_1 = function(model) {
  dd <- list(
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
  return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_0_2 = function(model) {
  up <- update_ratemat_indices

  ## TODO: spec problem -- numIters should have been kept in model
  numIters <- 3

  dd <- list(
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
  return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_0_4 = function(model) {
  dd <- list(
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
  return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_0_5 = function(model) {
  dd <- list(
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
  return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_0_6 = function(model) {
  dd <- list(
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
  return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_1_0 = function(model) {
  unpack(sum_indices)
  dd <- list(
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
  return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_1_1 = function(model) {
  unpack(sum_indices)
  init_tv_mult = integer(0L)
  if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult
  dd <- list(
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
  return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_1_2 = function(model) {
      unpack(sum_indices)
      init_tv_mult = integer(0L)
      if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult
      dd <- list(
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
      return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_2_0 = function(model) {

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

  dd <- list(
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
      tv_type_id = null_to_int0(schedule$tv_type_id),
      #tv_type = null_to_int0(NULL),

      sumidx = null_to_int0(sumidx),
      sumcount = null_to_int0(unname(sumcount)),
      summandidx = null_to_int0(summandidx),

      factr_spi = null_to_int0(factr_indices$spi_factr),
      factr_count = null_to_int0(factr_indices$count),
      factr_spi_compute = null_to_int0(factr_indices$spi),
      factr_modifier = null_to_int0(factr_indices$modifier),

      powidx = null_to_int0(pow_indices$powidx),
      powarg1idx = null_to_int0(pow_indices$powarg1idx),
      powarg2idx = null_to_int0(pow_indices$powarg2idx),
      powconstidx = null_to_int0(pow_indices$powconstidx),

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

      obs_do_sim_constraint = isTRUE(do_sim_constraint),
      obs_sim_lower_bound = null_to_num0(sim_lower_bound),

      numIterations = int0_to_0(null_to_0(iters))
    ),
    parameters = list(params = c(unlist(params)),
                      tv_mult = null_to_num0(init_tv_mult)),
    map = map,
    DLL = DLL,
    silent = getOption("MP_silent_tmb_function")
  )
  return(dd)
}

tmb_fun_args_by_spec.spec_ver_0_2_1 = function(model) {
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

  dd <- list(
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
      tv_type_id = null_to_int0(schedule$tv_type_id),
      #tv_type = null_to_int0(NULL),

      sumidx = null_to_int0(sumidx),
      sumcount = null_to_int0(unname(sumcount)),
      summandidx = null_to_int0(summandidx),

      factr_spi = null_to_int0(factr_indices$spi_factr),
      factr_count = null_to_int0(factr_indices$count),
      factr_spi_compute = null_to_int0(factr_indices$spi),
      factr_modifier = null_to_int0(factr_indices$modifier),

      powidx = null_to_int0(pow_indices$powidx),
      powarg1idx = null_to_int0(pow_indices$powarg1idx),
      powarg2idx = null_to_int0(pow_indices$powarg2idx),
      powconstidx = null_to_int0(pow_indices$powconstidx),

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
      lag_diff_delay_n = null_to_int0_mat(lag_diff$delay_n),

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

      obs_do_sim_constraint = isTRUE(do_sim_constraint),
      obs_sim_lower_bound = null_to_num0(sim_lower_bound),

      numIterations = int0_to_0(null_to_0(iters))
    ),
    parameters = list(params = c(unlist(params)),
                      tv_mult = null_to_num0(init_tv_mult)),
    map = map,
    DLL = DLL,
    silent = getOption("MP_silent_tmb_function")
  )
  return(dd)
}

#
# tmb_fun_args_by_spec.default <- function(model) {
#     ## make ad functions for different spec versions
#     if (spec_ver_eq("0.0.1")) {
#         dd <- list(
#             data = list(
#                 state = c(state),
#                 ratemat = ratemat,
#                 from = from,
#                 to = to,
#                 count = count,
#                 spi = spi,
#                 modifier = modifier
#             ),
#             parameters = list(params = c(unlist(params))),
#             DLL = DLL
#         )
#     } else if (spec_ver_eq("0.0.2")) {
#         up <- update_ratemat_indices
#
#         ## TODO: spec problem -- numIters should have been kept in model
#         numIters <- 3
#
#         dd <- list(
#             data = list(
#                 state = c(state),
#                 ratemat = ratemat,
#                 from = from,
#                 to = to,
#                 count = count,
#                 spi = spi,
#                 modifier = modifier,
#                 update_from = up$from,
#                 update_to = up$to,
#                 update_count = up$count,
#                 update_spi = up$spi,
#                 update_modifier = up$modifier,
#                 par_accum_indices = par_accum_indices,
#                 numIterations = numIters
#             ),
#             parameters = list(params = c(unlist(params))),
#             DLL = DLL
#         )
#     } else if (spec_ver_eq("0.0.4")) {
#         dd <- list(
#             data = list(
#                 state = c(state),
#                 ratemat = ratemat,
#                 from = from,
#                 to = to,
#                 count = count,
#                 spi = spi,
#                 modifier = modifier,
#                 updateidx = c(updateidx),
#                 breaks = breaks,
#                 count_of_tv_at_breaks = count_of_tv_at_breaks,
#                 tv_spi = schedule$tv_spi,
#                 tv_val = schedule$tv_val,
#                 par_accum_indices = par_accum_indices,
#                 numIterations = iters
#             ),
#             parameters = list(params = c(unlist(params))),
#             DLL = DLL
#         )
#     } else if (spec_ver_eq("0.0.5")) {
#         dd <- list(
#             data = list(
#                 state = c(state),
#                 ratemat = ratemat,
#                 from = from,
#                 to = to,
#                 count = count,
#                 spi = spi,
#                 modifier = modifier,
#                 updateidx = c(updateidx),
#                 breaks = breaks,
#                 count_of_tv_at_breaks = count_of_tv_at_breaks,
#                 tv_spi = schedule$tv_spi,
#                 tv_val = schedule$tv_val,
#                 par_accum_indices = par_accum_indices,
#                 do_hazard = do_hazard,
#                 numIterations = iters
#             ),
#             parameters = list(params = c(unlist(params))),
#             DLL = DLL
#         )
#     } else if (spec_ver_eq("0.0.6")) {
#         dd <- list(
#             data = list(
#                 state = c(state),
#                 ratemat = ratemat,
#                 from = from,
#                 to = to,
#                 count = count,
#                 spi = spi,
#                 modifier = modifier,
#                 updateidx = c(updateidx),
#                 breaks = breaks,
#                 count_of_tv_at_breaks = count_of_tv_at_breaks,
#                 tv_spi = schedule$tv_spi,
#                 tv_val = schedule$tv_val,
#                 tv_mult = schedule$Value,
#                 tv_orig = schedule$Type == "rel_orig",
#                 par_accum_indices = par_accum_indices,
#                 do_hazard = do_hazard,
#                 numIterations = iters
#             ),
#             parameters = list(params = c(unlist(params))),
#             DLL = DLL
#         )
#     } else if (spec_ver_eq("0.1.0")) {
#         unpack(sum_indices)
#         dd <- list(
#             data = list(
#                 state = c(state),
#                 ratemat = ratemat,
#                 from = from,
#                 to = to,
#                 count = count,
#                 spi = spi,
#                 modifier = modifier,
#                 updateidx = c(updateidx),
#                 breaks = breaks,
#                 count_of_tv_at_breaks = count_of_tv_at_breaks,
#                 tv_spi = schedule$tv_spi,
#                 tv_val = schedule$tv_val,
#                 tv_mult = schedule$Value,
#                 tv_orig = schedule$Type == "rel_orig",
#                 sumidx = sumidx,
#                 sumcount = unname(sumcount),
#                 summandidx = summandidx,
#                 par_accum_indices = par_accum_indices,
#                 do_hazard = do_hazard,
#                 numIterations = iters
#             ),
#             parameters = list(params = c(unlist(params))),
#             DLL = DLL
#         )
#     } else if (spec_ver_eq("0.1.1")) {
#         unpack(sum_indices)
#         init_tv_mult = integer(0L)
#         if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult
#         dd <- list(
#             data = list(
#                 state = c(state),
#                 ratemat = ratemat,
#                 from = null_to_int0(from),
#                 to = null_to_int0(to),
#                 count = null_to_int0(count),
#                 spi = null_to_int0(spi),
#                 modifier = null_to_int0(modifier),
#                 updateidx = null_to_int0(c(updateidx)),
#                 breaks = null_to_int0(breaks),
#                 count_of_tv_at_breaks = null_to_int0(count_of_tv_at_breaks),
#                 tv_val = null_to_num0(schedule$tv_val),
#                 tv_spi = null_to_int0(schedule$tv_spi),
#                 tv_spi_unique = null_to_int0(sort(unique(schedule$tv_spi))),
#                 # tv_mult = schedule$Value,  # moved to parameter vector
#                 tv_orig = null_to_log0(schedule$Type == "rel_orig"),
#                 tv_abs = null_to_log0(schedule$Type == "abs"),
#
#                 sumidx = null_to_int0(sumidx),
#                 sumcount = null_to_int0(unname(sumcount)),
#                 summandidx = null_to_int0(summandidx),
#
#                 do_make_state = isTRUE(do_make_state),
#                 max_iters_eig_pow_meth = int0_to_0(null_to_0(max_iters_eig_pow_meth)),
#                 tol_eig_pow_meth = null_to_num0(tol_eig_pow_meth),
#
#                 outflow_row_count = null_to_int0(outflow$row_count),
#                 outflow_col_count = null_to_int0(outflow$col_count),
#                 outflow_rows = null_to_int0(outflow$rows),
#                 outflow_cols = null_to_int0(outflow$cols),
#
#                 linearized_outflow_row_count = null_to_int0(linearized_outflow$row_count),
#                 linearized_outflow_col_count = null_to_int0(linearized_outflow$col_count),
#                 linearized_outflow_rows = null_to_int0(linearized_outflow$rows),
#                 linearized_outflow_cols = null_to_int0(linearized_outflow$cols),
#
#                 lin_param_vals = null_to_num0(linearized_params$lin_param_vals),
#                 lin_param_count = null_to_int0(linearized_params$lin_param_count),
#                 lin_param_idx = null_to_int0(linearized_params$lin_param_idx),
#
#                 df_state_par_idx = null_to_int0(disease_free$df_state_par_idx),
#                 df_state_count = null_to_int0(disease_free$df_state_count),
#                 df_state_idx = null_to_int0(disease_free$df_state_idx),
#
#                 im_all_drop_eigen_idx = null_to_int0(initialization_mapping$all_drop_eigen_idx),
#                 im_eigen_drop_infected_idx = null_to_int0(initialization_mapping$eigen_drop_infected_idx),
#                 im_all_to_infected_idx = null_to_int0(initialization_mapping$all_to_infected_idx),
#                 im_susceptible_idx = null_to_int0(initialization_mapping$susceptible_idx),
#
#                 ip_total_idx = int0_to_0(null_to_0(initial_population$total_idx)),
#                 ip_infected_idx = int0_to_0(null_to_0(initial_population$infected_idx)),
#
#                 do_hazard = isTRUE(do_hazard),
#                 do_hazard_lin = isTRUE(do_hazard_lin),
#                 do_approx_hazard = isTRUE(do_approx_hazard),
#                 do_approx_hazard_lin = isTRUE(do_approx_hazard_lin),
#                 haz_eps = haz_eps,
#
#                 numIterations = int0_to_0(null_to_0(iters))
#             ),
#             parameters = list(params = c(unlist(params)),
#                               tv_mult = init_tv_mult),
#             DLL = DLL
#         )
#     } else if (spec_ver_eq("0.1.2")) {
#       unpack(sum_indices)
#       init_tv_mult = integer(0L)
#       if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult
#       dd <- list(
#         data = list(
#           state = c(state),
#           ratemat = ratemat,
#           from = null_to_int0(from),
#           to = null_to_int0(to),
#           count = null_to_int0(count),
#           spi = null_to_int0(spi),
#           modifier = null_to_int0(modifier),
#           updateidx = null_to_int0(c(updateidx)),
#           breaks = null_to_int0(breaks),
#           count_of_tv_at_breaks = null_to_int0(count_of_tv_at_breaks),
#           tv_val = null_to_num0(schedule$tv_val),
#           tv_spi = null_to_int0(schedule$tv_spi),
#           tv_spi_unique = null_to_int0(sort(unique(schedule$tv_spi))),
#           # tv_mult = schedule$Value,  # moved to parameter vector
#           tv_orig = null_to_log0(schedule$Type == "rel_orig"),
#           tv_abs = null_to_log0(schedule$Type == "abs"),
#
#           sumidx = null_to_int0(sumidx),
#           sumcount = null_to_int0(unname(sumcount)),
#           summandidx = null_to_int0(summandidx),
#
#           factr_spi = null_to_int0(factr_indices$spi_factr),
#           factr_count = null_to_int0(factr_indices$count),
#           factr_spi_compute = null_to_int0(factr_indices$spi),
#           factr_modifier = null_to_int0(factr_indices$modifier),
#
#           do_make_state = isTRUE(do_make_state),
#           max_iters_eig_pow_meth = int0_to_0(null_to_0(max_iters_eig_pow_meth)),
#           tol_eig_pow_meth = null_to_num0(tol_eig_pow_meth),
#
#           outflow_row_count = null_to_int0(outflow$row_count),
#           outflow_col_count = null_to_int0(outflow$col_count),
#           outflow_rows = null_to_int0(outflow$rows),
#           outflow_cols = null_to_int0(outflow$cols),
#
#           linearized_outflow_row_count = null_to_int0(linearized_outflow$row_count),
#           linearized_outflow_col_count = null_to_int0(linearized_outflow$col_count),
#           linearized_outflow_rows = null_to_int0(linearized_outflow$rows),
#           linearized_outflow_cols = null_to_int0(linearized_outflow$cols),
#
#           lin_param_vals = null_to_num0(linearized_params$lin_param_vals),
#           lin_param_count = null_to_int0(linearized_params$lin_param_count),
#           lin_param_idx = null_to_int0(linearized_params$lin_param_idx),
#
#           df_state_par_idx = null_to_int0(disease_free$df_state_par_idx),
#           df_state_count = null_to_int0(disease_free$df_state_count),
#           df_state_idx = null_to_int0(disease_free$df_state_idx),
#
#           im_all_drop_eigen_idx = null_to_int0(initialization_mapping$all_drop_eigen_idx),
#           im_eigen_drop_infected_idx = null_to_int0(initialization_mapping$eigen_drop_infected_idx),
#           im_all_to_infected_idx = null_to_int0(initialization_mapping$all_to_infected_idx),
#           im_susceptible_idx = null_to_int0(initialization_mapping$susceptible_idx),
#
#           ip_total_idx = int0_to_0(null_to_0(initial_population$total_idx)),
#           ip_infected_idx = int0_to_0(null_to_0(initial_population$infected_idx)),
#
#           do_hazard = isTRUE(do_hazard),
#           do_hazard_lin = isTRUE(do_hazard_lin),
#           do_approx_hazard = isTRUE(do_approx_hazard),
#           do_approx_hazard_lin = isTRUE(do_approx_hazard_lin),
#           # haz_eps = haz_eps,
#
#           sri_output = null_to_int0(sim_report_expr_indices$sri_output),
#           sr_count = null_to_int0(sim_report_expr_indices$sr_count),
#           sri = null_to_int0(sim_report_expr_indices$sri),
#           sr_modifier = null_to_int0(sim_report_expr_indices$sr_modifier),
#
#           lag_diff_sri = null_to_int0(lag_diff$sri),
#           lag_diff_delay_n = null_to_int0(lag_diff$delay_n),
#
#           conv_sri = null_to_int0(conv$sri),
#           conv_c_prop_idx = null_to_int0(conv$c_prop_idx),
#           conv_c_delay_cv_idx = null_to_int0(conv$c_delay_cv_idx),
#           conv_c_delay_mean_idx = null_to_int0(conv$c_delay_mean_idx),
#           conv_qmax = null_to_int0(conv$qmax),
#
#           numIterations = int0_to_0(null_to_0(iters))
#         ),
#         parameters = list(params = c(unlist(params)),
#                           tv_mult = init_tv_mult),
#         DLL = DLL
#       )
#
#     } else if (spec_ver_eq("0.2.0")) {
#
#       unpack(sum_indices)
#       unpack(opt_params)
#
#       # update parameter vectors -----------------
#       init_tv_mult = integer(0L)
#       if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult
#
#       if(isTRUE(exists_opt_params(model))) {
#
#         # tell tmb what parameters to put in the objective function
#         map = ad_fun_map
#
#         # pass transformed parameters to the objective function
#         params = tmb_params_trans(model, vec_type = 'params')
#         init_tv_mult = tmb_params_trans(model, vec_type = 'tv_mult')
#
#       } else {
#         map = list()
#       }
#       # -------------------------------------------
#
#       if (!all(between(observed$time_step, 2, iters + 1))) {
#         stop('observations outside of simulation range')
#       }
#
#       dd <- list(
#         data = list(
#           state = c(state),
#           ratemat = ratemat,
#           from = null_to_int0(from),
#           to = null_to_int0(to),
#           count = null_to_int0(count),
#           spi = null_to_int0(spi),
#           modifier = null_to_int0(modifier),
#           updateidx = null_to_int0(c(updateidx)),
#           breaks = null_to_int0(breaks),
#           count_of_tv_at_breaks = null_to_int0(count_of_tv_at_breaks),
#
#           tv_val = null_to_num0(schedule$tv_val),
#           tv_spi = null_to_int0(schedule$tv_spi),
#           tv_spi_unique = null_to_int0(sort(unique(schedule$tv_spi))),
#           # tv_mult = schedule$Value,  # moved to parameter vector
#           tv_orig = null_to_log0(schedule$Type == "rel_orig"),
#           tv_abs = null_to_log0(schedule$Type == "abs"),
#           tv_type_id = null_to_int0(schedule$tv_type_id),
#           #tv_type = null_to_int0(NULL),
#
#           sumidx = null_to_int0(sumidx),
#           sumcount = null_to_int0(unname(sumcount)),
#           summandidx = null_to_int0(summandidx),
#
#           factr_spi = null_to_int0(factr_indices$spi_factr),
#           factr_count = null_to_int0(factr_indices$count),
#           factr_spi_compute = null_to_int0(factr_indices$spi),
#           factr_modifier = null_to_int0(factr_indices$modifier),
#
#           powidx = null_to_int0(pow_indices$powidx),
#           powarg1idx = null_to_int0(pow_indices$powarg1idx),
#           powarg2idx = null_to_int0(pow_indices$powarg2idx),
#           powconstidx = null_to_int0(pow_indices$powconstidx),
#
#           do_make_state = isTRUE(do_make_state),
#           max_iters_eig_pow_meth = int0_to_0(null_to_0(max_iters_eig_pow_meth)),
#           tol_eig_pow_meth = null_to_num0(tol_eig_pow_meth),
#
#           outflow_row_count = null_to_int0(outflow$row_count),
#           outflow_col_count = null_to_int0(outflow$col_count),
#           outflow_rows = null_to_int0(outflow$rows),
#           outflow_cols = null_to_int0(outflow$cols),
#
#           linearized_outflow_row_count = null_to_int0(linearized_outflow$row_count),
#           linearized_outflow_col_count = null_to_int0(linearized_outflow$col_count),
#           linearized_outflow_rows = null_to_int0(linearized_outflow$rows),
#           linearized_outflow_cols = null_to_int0(linearized_outflow$cols),
#
#           lin_param_vals = null_to_num0(linearized_params$lin_param_vals),
#           lin_param_count = null_to_int0(linearized_params$lin_param_count),
#           lin_param_idx = null_to_int0(linearized_params$lin_param_idx),
#
#           df_state_par_idx = null_to_int0(disease_free$df_state_par_idx),
#           df_state_count = null_to_int0(disease_free$df_state_count),
#           df_state_idx = null_to_int0(disease_free$df_state_idx),
#
#           im_all_drop_eigen_idx = null_to_int0(initialization_mapping$all_drop_eigen_idx),
#           im_eigen_drop_infected_idx = null_to_int0(initialization_mapping$eigen_drop_infected_idx),
#           im_all_to_infected_idx = null_to_int0(initialization_mapping$all_to_infected_idx),
#           im_susceptible_idx = null_to_int0(initialization_mapping$susceptible_idx),
#
#           ip_total_idx = int0_to_0(null_to_0(initial_population$total_idx)),
#           ip_infected_idx = int0_to_0(null_to_0(initial_population$infected_idx)),
#
#           do_hazard = isTRUE(do_hazard),
#           do_hazard_lin = isTRUE(do_hazard_lin),
#           do_approx_hazard = isTRUE(do_approx_hazard),
#           do_approx_hazard_lin = isTRUE(do_approx_hazard_lin),
#           # haz_eps = haz_eps,
#
#           sri_output = null_to_int0(sim_report_expr_indices$sri_output),
#           sr_count = null_to_int0(sim_report_expr_indices$sr_count),
#           sri = null_to_int0(sim_report_expr_indices$sri),
#           sr_modifier = null_to_int0(sim_report_expr_indices$sr_modifier),
#
#           lag_diff_sri = null_to_int0(lag_diff$sri),
#           lag_diff_delay_n = null_to_int0(lag_diff$delay_n),
#
#           conv_sri = null_to_int0(conv$sri),
#           conv_c_prop_idx = null_to_int0(conv$c_prop_idx),
#           conv_c_delay_cv_idx = null_to_int0(conv$c_delay_cv_idx),
#           conv_c_delay_mean_idx = null_to_int0(conv$c_delay_mean_idx),
#           conv_qmax = null_to_int0(conv$qmax),
#
#           obs_var_id = null_to_int0(observed$variable_id),
#           obs_loss_id = null_to_int0(observed$loss_id),
#           obs_loss_param_count = null_to_int0(observed$loss_param_count),
#           obs_spi_loss_param = null_to_int0(observed$spi_loss_param),
#           obs_time_step = null_to_int0(observed$time_step),
#           obs_history_col_id = null_to_int0(observed$history_col_id),
#           obs_value = observed$observed, # don't need to worry about missing values because they are omitted
#
#
#           opt_param_id = null_to_int0(index_table$opt_param_id),
#           opt_trans_id = null_to_int0(index_table$param_trans_id),
#           opt_count_reg_params = null_to_int0(index_table$count_hyperparams),
#           opt_reg_params = null_to_num0(hyperparameters),
#           opt_reg_family_id = null_to_int0(index_table$prior_distr_id),
#
#           opt_tv_param_id = null_to_int0(index_tv_table$opt_tv_mult_id),
#           opt_tv_trans_id = null_to_int0(index_tv_table$param_trans_id),
#           opt_tv_count_reg_params = null_to_int0(index_tv_table$count_hyperparams),
#           opt_tv_reg_params = null_to_num0(hyperparameters_tv),
#           opt_tv_reg_family_id = null_to_int0(index_tv_table$prior_distr_id),
#
#           obs_do_sim_constraint = isTRUE(do_sim_constraint),
#           obs_sim_lower_bound = null_to_num0(sim_lower_bound),
#
#           numIterations = int0_to_0(null_to_0(iters))
#         ),
#         parameters = list(params = c(unlist(params)),
#                           tv_mult = null_to_num0(init_tv_mult)),
#         map = map,
#         DLL = DLL,
#         silent = getOption("MP_silent_tmb_function")
#       )
#     } else if (spec_ver_eq('0.2.1')) {
#
#       unpack(sum_indices)
#       unpack(opt_params)
#
#       # update parameter vectors -----------------
#       init_tv_mult = integer(0L)
#       if(!is.null(schedule$init_tv_mult)) init_tv_mult = schedule$init_tv_mult
#
#       if(isTRUE(exists_opt_params(model))) {
#
#         # tell tmb what parameters to put in the objective function
#         map = ad_fun_map
#
#         # pass transformed parameters to the objective function
#         params = tmb_params_trans(model, vec_type = 'params')
#         init_tv_mult = tmb_params_trans(model, vec_type = 'tv_mult')
#
#       } else {
#         map = list()
#       }
#       # -------------------------------------------
#
#       if (!all(between(observed$time_step, 2, iters + 1))) {
#         stop('observations outside of simulation range')
#       }
#
#       dd <- list(
#         data = list(
#           state = c(state),
#           ratemat = ratemat,
#           from = null_to_int0(from),
#           to = null_to_int0(to),
#           count = null_to_int0(count),
#           spi = null_to_int0(spi),
#           modifier = null_to_int0(modifier),
#           updateidx = null_to_int0(c(updateidx)),
#           breaks = null_to_int0(breaks),
#           count_of_tv_at_breaks = null_to_int0(count_of_tv_at_breaks),
#
#           tv_val = null_to_num0(schedule$tv_val),
#           tv_spi = null_to_int0(schedule$tv_spi),
#           tv_spi_unique = null_to_int0(sort(unique(schedule$tv_spi))),
#           # tv_mult = schedule$Value,  # moved to parameter vector
#           tv_orig = null_to_log0(schedule$Type == "rel_orig"),
#           tv_abs = null_to_log0(schedule$Type == "abs"),
#           tv_type_id = null_to_int0(schedule$tv_type_id),
#           #tv_type = null_to_int0(NULL),
#
#           sumidx = null_to_int0(sumidx),
#           sumcount = null_to_int0(unname(sumcount)),
#           summandidx = null_to_int0(summandidx),
#
#           factr_spi = null_to_int0(factr_indices$spi_factr),
#           factr_count = null_to_int0(factr_indices$count),
#           factr_spi_compute = null_to_int0(factr_indices$spi),
#           factr_modifier = null_to_int0(factr_indices$modifier),
#
#           powidx = null_to_int0(pow_indices$powidx),
#           powarg1idx = null_to_int0(pow_indices$powarg1idx),
#           powarg2idx = null_to_int0(pow_indices$powarg2idx),
#           powconstidx = null_to_int0(pow_indices$powconstidx),
#
#           do_make_state = isTRUE(do_make_state),
#           max_iters_eig_pow_meth = int0_to_0(null_to_0(max_iters_eig_pow_meth)),
#           tol_eig_pow_meth = null_to_num0(tol_eig_pow_meth),
#
#           outflow_row_count = null_to_int0(outflow$row_count),
#           outflow_col_count = null_to_int0(outflow$col_count),
#           outflow_rows = null_to_int0(outflow$rows),
#           outflow_cols = null_to_int0(outflow$cols),
#
#           linearized_outflow_row_count = null_to_int0(linearized_outflow$row_count),
#           linearized_outflow_col_count = null_to_int0(linearized_outflow$col_count),
#           linearized_outflow_rows = null_to_int0(linearized_outflow$rows),
#           linearized_outflow_cols = null_to_int0(linearized_outflow$cols),
#
#           lin_param_vals = null_to_num0(linearized_params$lin_param_vals),
#           lin_param_count = null_to_int0(linearized_params$lin_param_count),
#           lin_param_idx = null_to_int0(linearized_params$lin_param_idx),
#
#           df_state_par_idx = null_to_int0(disease_free$df_state_par_idx),
#           df_state_count = null_to_int0(disease_free$df_state_count),
#           df_state_idx = null_to_int0(disease_free$df_state_idx),
#
#           im_all_drop_eigen_idx = null_to_int0(initialization_mapping$all_drop_eigen_idx),
#           im_eigen_drop_infected_idx = null_to_int0(initialization_mapping$eigen_drop_infected_idx),
#           im_all_to_infected_idx = null_to_int0(initialization_mapping$all_to_infected_idx),
#           im_susceptible_idx = null_to_int0(initialization_mapping$susceptible_idx),
#
#           ip_total_idx = int0_to_0(null_to_0(initial_population$total_idx)),
#           ip_infected_idx = int0_to_0(null_to_0(initial_population$infected_idx)),
#
#           do_hazard = isTRUE(do_hazard),
#           do_hazard_lin = isTRUE(do_hazard_lin),
#           do_approx_hazard = isTRUE(do_approx_hazard),
#           do_approx_hazard_lin = isTRUE(do_approx_hazard_lin),
#           # haz_eps = haz_eps,
#
#           sri_output = null_to_int0(sim_report_expr_indices$sri_output),
#           sr_count = null_to_int0(sim_report_expr_indices$sr_count),
#           sri = null_to_int0(sim_report_expr_indices$sri),
#           sr_modifier = null_to_int0(sim_report_expr_indices$sr_modifier),
#
#           lag_diff_sri = null_to_int0(lag_diff$sri),
#           lag_diff_delay_n = null_to_int0_mat(lag_diff$delay_n),
#
#           conv_sri = null_to_int0(conv$sri),
#           conv_c_prop_idx = null_to_int0(conv$c_prop_idx),
#           conv_c_delay_cv_idx = null_to_int0(conv$c_delay_cv_idx),
#           conv_c_delay_mean_idx = null_to_int0(conv$c_delay_mean_idx),
#           conv_qmax = null_to_int0(conv$qmax),
#
#           obs_var_id = null_to_int0(observed$variable_id),
#           obs_loss_id = null_to_int0(observed$loss_id),
#           obs_loss_param_count = null_to_int0(observed$loss_param_count),
#           obs_spi_loss_param = null_to_int0(observed$spi_loss_param),
#           obs_time_step = null_to_int0(observed$time_step),
#           obs_history_col_id = null_to_int0(observed$history_col_id),
#           obs_value = observed$observed, # don't need to worry about missing values because they are omitted
#
#
#           opt_param_id = null_to_int0(index_table$opt_param_id),
#           opt_trans_id = null_to_int0(index_table$param_trans_id),
#           opt_count_reg_params = null_to_int0(index_table$count_hyperparams),
#           opt_reg_params = null_to_num0(hyperparameters),
#           opt_reg_family_id = null_to_int0(index_table$prior_distr_id),
#
#           opt_tv_param_id = null_to_int0(index_tv_table$opt_tv_mult_id),
#           opt_tv_trans_id = null_to_int0(index_tv_table$param_trans_id),
#           opt_tv_count_reg_params = null_to_int0(index_tv_table$count_hyperparams),
#           opt_tv_reg_params = null_to_num0(hyperparameters_tv),
#           opt_tv_reg_family_id = null_to_int0(index_tv_table$prior_distr_id),
#
#           obs_do_sim_constraint = isTRUE(do_sim_constraint),
#           obs_sim_lower_bound = null_to_num0(sim_lower_bound),
#
#           numIterations = int0_to_0(null_to_0(iters))
#         ),
#         parameters = list(params = c(unlist(params)),
#                           tv_mult = null_to_num0(init_tv_mult)),
#         map = map,
#         DLL = DLL,
#         silent = getOption("MP_silent_tmb_function")
#       )
#     }
#     else {
#         stop("Construction of TMB functions is not supported by your installation of MacPan")
#     }
#     return(dd)
# }

##' Update Initial State
##'
##' Update the initial state of a \code{\link{flexmodel}} object
##' using an eigenvector-based approach to finding a near-disease-free
##' state
##'
##' @param model \code{\link{flexmodel}} object
##' @param silent warn if the model is not properly defined for
##' updating the initial state?
##'
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

# observed data and observation error ---------------------------------

#' @importFrom dplyr bind_rows
#' @rdname update_error_dist
#' @export
add_error_dist = function(model, ...) {
  new_loss_params = (list(...)
   %>% lapply(
     parse_and_resolve_loss_form,
     condensed_sim_report_names(model),
     names(model$params)
    )
   %>% bind_rows
  )
  new_loss_params = bind_rows(
    model$observed$loss_params,
    new_loss_params
  )
  update_loss_params(model, new_loss_params)
}

#' Update Error Distribution
#'
#' @param model \code{\link{flexmodel}} object
#' @param ... optional error distribution formulas. if no formulas are provided
#' the existing error distributions are removed.
#' @export
update_error_dist = function(model, ...) {
  if (length(list(...)) == 0L) {
    model$observed$loss_params = init_observed$loss_params
    return(model)
  }
  new_loss_params = (list(...)
   %>% lapply(
     parse_and_resolve_loss_form,
     condensed_sim_report_names(model),
     names(model$params)
    )
   %>% bind_rows
  )
  update_loss_params(model, new_loss_params)
}

#' @rdname update_error_dist
#' @export
reset_error_dist = function(model) {
  update_error_dist(model)
}

#' Update Observation Error
#'
#' @param model \code{\link{flexmodel}} object
#' @param loss_params TODO
#' @param regenerate_rates TODO
#' @export
update_loss_params = function(model, loss_params, regenerate_rates = TRUE) {

  dist_nms = unique(loss_params$Distribution)
  invalid_dist_nms = dist_nms[!dist_nms %in% valid_loss_functions]
  if (length(invalid_dist_nms) != 0L) {
    stop(
      'the following distributions were requested but not valid: ',
      paste0(invalid_dist_nms, collapse = ", ")
    )
  }

  model$observed$loss_params = loss_params

  if (regenerate_rates) {
    model = regen_model(model)
  }

  # HACK: refactor so that update_loss_params has a generic,
  #       so that this class setting becomes unconditional
  if (isTRUE(all.equal(model$observed$data, init_observed$data))) {
    class(model) = c("flexmodel_obs_error", "flexmodel")
  } else {
    class(model) = c(
      "flexmodel_to_calibrate",
      "flexmodel"
    )
  }
  model
}

#' Update Observed Data
#'
#' Attach a data set to a \code{\link{flexmodel}}
#' object. Any existing data will be removed.
#'
#' @param model \code{\link{flexmodel}} object
#' @param data observed data frame in long format to
#' compare with simulated trajectories. must have the following
#' columns: \code{date}, \code{var}, \code{value}.
#' @param loss_params TODO
#' @param regenerate_rates should rates be regenerated to sort out any
#' index misalignment that may be introduced by editing the model?
#'
#' @return a \code{\link{flexmodel}} object with an \code{observed} element,
#' which is a list with two elements: (1) \code{data} (the attached data)
#' and (2) \code{loss} (a data frame describing the additional parameters
#' that are required to fit the model to data -- these parameters are added
#' to the \code{params} of the model)
#'
#' @export
update_observed = function(model, data, loss_params = NULL, regenerate_rates = TRUE) {
  stopifnot(isTRUE(all.equal(c(names(data)), c("date", "var", "value"))))
  obsvars = unique(data$var)
  model$observed$data = data

  if (is.null(loss_params)) {
    if (nrow(model$observed$loss_params) == 0L) {
      # message(
      #   "\na default negative binomial error distribution for all observed\n",
      #   "variables has been assumed. one dispersion parameter for each\n",
      #   "variable has been added to the parameter vector and set to a\n",
      #   "default value of 1. it is recommended that update_error_dist is\n",
      #   "explicitly called."
      # )
      model$observed$loss_params = data.frame(
        Parameter = "nb_disp" %_% obsvars, # default choice: dispersion
        Distribution = "negative_binomial",   # default choice: negative binomial
        Variable = obsvars
      )
      model$params = expand_params_nb_disp(model$params, obsvars)
    } #else {
      # warning(
      #   "\nan error distribution was inherited from a modified flexmodel, \n",
      #   "and has not been explicitly modified to be consistent with the \n",
      #   "new observed data. it is recommended that update_error_dist is \n",
      #   "explicitly called."
      # )
    #}
  } else {
    model = update_loss_params(model, loss_params)
  }
  if (regenerate_rates) {
    model = regen_model(model)
  }
  class(model) = c(
    'flexmodel_to_calibrate',
    'flexmodel'
  )
  return(model)
}

# time variation updates ------------------

#' Update the Piece-Wise Parameter Time-Variation Schedule
#'
#' @param model \code{\link{flexmodel}} object
#' @param params_timevar data frame with scheduling for piece-wise
#' constant parameter variation (TODO: direct to other help pages)
#' @param regenerate_rates should the rates in the flexmodel be
#' regenerated? If you don't know what you are doing, you should
#' select \code{TRUE}, which is the default.
#' @export
update_piece_wise = function(model, params_timevar, regenerate_rates = TRUE) {
  spec_check(
    introduced_version = '0.0.4',
    feature = 'piece-wise time variation of parameters'
  )

  #
  schedule <- (params_timevar
   %>% mutate(Date = as.Date(Date))
   %>% mutate(breaks = (Date
      %>% difftime(model$start_date, units = 'days')
      %>% as.integer
   ))
   %>% mutate(tv_spi = find_vec_indices(Symbol, c(model$state, model$params)))
   %>% mutate(tv_type_id = as.numeric(factor(Type, levels = valid_tv_types)))
   %>% arrange(breaks, tv_spi)
  )

    if(spec_ver_gt("0.1.0")) {
      schedule = (schedule
        %>% mutate(init_tv_mult = replace(
          Value,
          which(is.na(Value)),
          1
        ))
        %>% mutate(last_tv_mult = init_tv_mult)
      )
    }

    count_of_tv_at_breaks <- c(table(schedule$breaks))

    if (nrow(schedule) > 0L) {
      schedule$tv_val <- NA
    }
    new_param <- TRUE
    ns <- nrow(schedule)
    if (!all(schedule$Type %in% valid_tv_types)) {
      stop(
        "Unrecognized break-point type.\n",
        "Only the following are allowed:\n",
        paste0(valid_tv_types, collapse = ', ')
      )
    }
    for (i in seq_len(ns)) {
      if ((new_param | grepl("^rel_orig", schedule$Type[i])) & (schedule$Type[i] != "abs")) {
        old_val <- model$params[schedule$Symbol[i]]
      } else if (grepl("^rel_prev", schedule$Type[i])) {
        old_val <- schedule$tv_val[i - 1]
      } else if (schedule$Type[i] == "abs") {
        old_val <- 1
      } else {
        stop(
          "parameter time-variation functionality is broken. ",
          "please contact package maintainers."
        )
      }
      if (grepl("_logit$", schedule$Type[i])) {
        schedule$tv_val[i] <- plogis(
          qlogis(old_val) + schedule$Value[i]
        )
      } else {
        schedule$tv_val[i] <- old_val * schedule$Value[i]
      }
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
    if (regenerate_rates) {
      model = regen_model(model)
    }
    if (spec_ver_gt('0.1.2')) {
      model$timevar$piece_wise$schedule = (model$timevar$piece_wise$schedule
       %>% mutate(tv_mult_id = seq_len(nrow(model$timevar$piece_wise$schedule)))
      )
    }
    model
}

#' @rdname update_piece_wise
#' @export
add_piece_wise = function(model, params_timevar, regenerate_rates = TRUE) {
  tv_cols = c("Date", "Symbol", "Value", "Type")
  if (!is.null(model$model_to_calibrate)) {
    ptv_to_cal = rbind(
      model$model_to_calibrate$timevar$piece_wise$schedule[tv_cols],
      params_timevar
    )
    model$model_to_calibrate = update_piece_wise(
      model$model_to_calibrate,
      ptv_to_cal,
      regenerate_rates
    )
  }
  ptv_to_cal = rbind(
    model$timevar$piece_wise$schedule[tv_cols],
    params_timevar
  )
  update_piece_wise(
    model,
    ptv_to_cal,
    regenerate_rates
  )
}

#' @rdname update_piece_wise
#' @export
initialize_piece_wise = function(model) {
  model$opt_tv_params = NULL
  model$timevar = list(
    piece_wise = list(
      breaks = integer(0L),
      count_of_tv_at_breaks = integer(0L),
      schedule = data.frame(
        Date = as.Date(numeric(0L), origin = "1970-01-01"),
        Symbol = character(0L),
        Value = numeric(0L),
        Type = character(0L),
        breaks = integer(0L),
        tv_spi = integer(0L),
        tv_type_id = integer(0L),
        init_tv_mult = numeric(0L),
        last_tv_mult = numeric(0L),
        tv_val = numeric(0L),
        tv_mult_id = integer(0L)
      )
    )
  )
  model
}


# param updates -----------------

#' Update Default Parameters
#'
#' @param model \code{\link{flexmodel}} object
#' @param ... named vectors of parameter names
#' @seealso \code{\link{pars_base_sim<-}}
#' @export
update_params = function(model, ...) {
  params_update = unlist(list(...))
  stopifnot(!is.null(names(params_update)))
  stopifnot(!any(duplicated(names(params_update))))
  stopifnot(all(grepl(
    wrap_exact(getOption("MP_name_search_regex")),
    names(params_update)
  )))
  model$params[names(params_update)] = params_update
  regen_model(model)
}

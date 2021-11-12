##' Initialize Flexible Compartmental Model
##'
##' https://canmod.net/misc/flex_specs
##'
##' @param params a \code{param_pansim} object
##' @param state a \code{state_pansim} object
##' @param struc a \code{struc_pansim} object
##' @param start_date simulation start date
##' @param end_date simulation end date
##' @param params_timevar data frame with scheduling for piece-wise
##' constant parameter variation (TODO: direct to other help pages)
##' @param do_hazard should hazard simulation steps be used?
##' (https://canmod.net/misc/flex_specs#v0.0.5) -- only used
##' if \code{spec_ver_gt('0.0.4')}
##' @family flexmodels
##' @return flexmodel object representing a compartmental model
##' @export
init_model <- function(params, state, struc = NULL,
                       start_date = NULL, end_date = NULL,
                       params_timevar = NULL,
                       do_hazard = TRUE, ...) {
    check_spec_ver_archived()
    name_regex = paste0("^", getOption("MP_name_search_regex"), "$")
    if(!all(grepl(name_regex, c(names(params), names(state))))) {
        stop("only syntactically valid r names can be used ",
             "for state variables and parameters")
    }

    # need to do this before we ruin the pansim structure
    ratemat = make_ratemat(state, params, sparse = TRUE)

    # TODO: keep an eye on this -- i think that the
    # flex framework should _not_ use parameter
    # lists and rather stick to numeric vectors, but
    # not totally sure
    params = setNames(unlist(params), names(params))

    model <- list(
        state = state,
        params = params,
        ratemat = ratemat,
        rates = list(),
        name_regex = name_regex
    )

    if (spec_ver_gt("0.0.1")) {
        model$parallel_accumulators <- character(0L)
    }

    if ((!is.null(start_date)) & (!is.null(end_date))) {
        spec_check(introduced_version = "0.0.3",
                   feature = "Start and end dates")
        model$start_date <- as.Date(start_date)
        model$end_date <- as.Date(end_date)
        model$iters <- as.integer(model$end_date - model$start_date)
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
            schedule <- (
                params_timevar
                    %>% mutate(Date = as.Date(Date))
                    %>% mutate(breaks = as.integer(Date - model$start_date))
                    %>% mutate(tv_spi = find_vec_indices(
                        Symbol, c(state, params)))
                    %>% arrange(breaks, tv_spi)
            )
            count_of_tv_at_breaks <- c(table(schedule$breaks))

            schedule$tv_val <- NA
            new_param <- TRUE
            ns <- nrow(schedule)
            for (i in 1:ns) {
                if (new_param | schedule$Type[i] == "rel_orig") {
                    old_val <- params[schedule$Symbol[i]]
                } else if (schedule$Type[i] == "rel_prev") {
                    old_val <- schedule$tv_val[i - 1]
                } else {
                    stop(
                        "Unrecognized break-point type.\n",
                        "Only rel_orig and rel_prev are allowed"
                    )
                }
                schedule$tv_val[i] <- old_val * schedule$Value[i]
                new_param <- schedule$Symbol[i] != schedule$Symbol[min(ns, i + i)]
            }

            attr(model$params, "tv_param_indices") <- setNames(find_vec_indices(
                unique(schedule$Symbol),
                params
            ), unique(schedule$Symbol))


            model$timevar$piece_wise <- list(
                ## schedule includes tv_spi and tv_val as in spec
                schedule = schedule,
                breaks = as.integer(names(count_of_tv_at_breaks)),
                count_of_tv_at_breaks = unname(count_of_tv_at_breaks)
            )
        } ## >v0.0.3
    } else {
        model$timevar = list(
            piece_wise = list(
                breaks = integer(0L),
                count_of_tv_at_breaks = integer(0L),
                schedule = list(
                    tv_spi = integer(0L),
                    tv_val = numeric(0L),
                    Value = numeric(0L),
                    Type = character(0L)
                )
            )
        )
    }
    if (spec_ver_gt("0.0.4")) model$do_hazard <- do_hazard
    if (spec_ver_gt("0.0.6")) {
        model$sums = list()
        model$sum_vector = numeric(0L)
    }

    if (spec_ver_eq("0.1.1")) {
        model$disease_free = list(
            state = list(
                simple = list(),
                state_mappings = list()
            )
        )
        model$linearized_params = list()
        model$outflow = list()
        model$linearized_outflow = list()
    }

    model$tmb_indices <- list(
        make_ratemat_indices = list(
            from = integer(0L),
            to = integer(0L),
            count = integer(0L),
            spi = integer(0L),
            modifier = integer(0L)
        ),
        par_accum_indices = integer(0L),
        updateidx = integer(0L),
        sum_indices = list(
            sumidx = integer(0L),
            sumcount = integer(0L),
            summandidx = integer(0L)
        )
    )

    structure(model, class = "flexmodel")
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
##' @family flexmodels
##' @export
rate <- function(from, to, formula, state, params, sums, ratemat) {
    ## TODO: test for formula structure
    ## TODO: test that from and to are available in the state vector
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

    product_list <- function(x) {
        x$factors <- (x$formula
                      %>% parse_formula
                      %>% lapply(factor_table) %>% bind_rows(.id = "prod_indx")
        )
        x$ratemat_indices <-
            do.call(McMasterPandemic:::pfun, c(x[c("from", "to")], list(mat = M)))
        if(nrow(x$ratemat_indices) > 1L) {
            stop('you are referring to more than one element of the rate matrix\n',
                 'try using rep_rate instead of rate')
        }
        x$factors$var_indx <- find_vec_indices(
            x$factors$var,
            c(state, params, sums))

        missing_vars = x$factors$var[sapply(x$factors$var_indx, length) == 0L]
        if(length(missing_vars) > 0L) {
            stop("The following variables were used to define the model,\n",
                 "but they could not be found in the state, parameter or sum vectors:\n",
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
    added_rate <- (
        rate(from, to, formula, model$state, model$params, model$sum_vector, model$ratemat)
            %>% list()
            %>% setNames(paste(from, to, sep = "_to_"))
    )
    model$rates <- c(model$rates, added_rate)
    return(model)
}

#' Repeat a Rate for Several Rate Matrix Elements
#'
#' @param model flexmodel
#' @param indices two-column matrix of indices with column names "from_pos"
#' and "to_pos", locaing elements of the rate matrix in model
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
        MoreArgs = nlist(formula, state, params, sums, ratemat),
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

    map_fun = switch(
        match.arg(mapping),
        pairwise = pwise,
        blockwise = block)

    unpack(model)

    indices = map_fun(from, to, ratemat)


    from = rownames(ratemat)[indices[,'from_pos']]
    to = colnames(ratemat)[indices[,'to_pos']]

    lst = mapply(rate, from, to, as.character(formula),
                 MoreArgs = nlist(state, params, sums, ratemat),
                 SIMPLIFY = FALSE, USE.NAMES = FALSE)
    nms = mapply(paste, from, to, MoreArgs = list(sep = "_to_"))
    model$rates <- c(rates, setNames(lst, nms))

    return(model)
}

#' Specify Matrix of Rates
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

# sums of state variables and parameters ---------------------

#' @param sum name of sum of state variables and parameters
#' @param summands character vector of regular expressions for identifying
#' state variables and parameters to sum together
#' @param state pansim_state object
#' @param params pansim_param object
#' @return indices of summands
#' @family flexmodels
#' @export
state_param_sum = function(sum, summands, state, params) {
    spec_check('0.1.0', 'sums of state variables and parameters')
    sp = c(state, params)
    summands = (summands
        %>% lapply(grep, names(sp), value = TRUE)
        %>% unlist
    )
    ii = find_vec_indices(summands, sp)
    val = sum(sp[ii])
    list(
        summands = summands,
        sum_indices = ii,
        initial_value = val)
}

#' @rdname state_param_sum
#' @inheritParams state_param_sum
#' @param model flexmodel
#' @family flexmodels
#' @export
add_state_param_sum = function(model, sum, summands) {
    model$sums[[sum]] = state_param_sum(
        sum, summands, model$state, model$params)

    # assumes that order of sums doesn't change!
    model$sum_vector = get_sum_initial_value(model) %>% unlist
    model
}

# parallel accumulators ---------------------

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
    (state_patterns
        %>% lapply(function(x) {
            grep(x, colnames(model$ratemat), value = TRUE)
        })
        %>% unlist()
    )
}

##' @family flexmodels
##' @export
add_linearized_outflow = function(model, state_patterns, flow_state_patterns) {
    model$linearized_outflow = append(
        model$linearized_outflow,
        list(outflow(model, state_patterns, flow_state_patterns)))
    return(model)
}

##' @family flexmodels
##' @export
add_outflow = function(model, state_patterns, flow_state_patterns) {
    model$outflow = append(
        model$outflow,
        list(outflow(model, state_patterns, flow_state_patterns)))
    return(model)
}

##' @family flexmodels
##' @export
outflow = function(model, state_patterns, flow_state_patterns) {
    nlist(state_patterns, flow_state_patterns)
}

##' @family
##' @export
add_eigen_scaler = function(model, state_name) {
    model$eigen_scaler = state_name
    model
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
    model$disease_free$state$simple = c(
        model$disease_free$state$simple,
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
##' \code{infected} is the scaler of the normalized infected components of the
##' eigenvector. \code{total - infected} is multiplied by
##' \code{1/length(initial_susceptible)}
##'
##' In the future we might want the ability to specify a distribution
##' for the susceptible distribution, but not currently.
##'
##' @param total name of a single parameter to represent the total size of the
##' population -- over all compartments
##' @param infected name of a single parameter to represent the total size of
##' the infected population -- over all infected compartments
##' @family flexmodels
##' @export
initial_population = function(model, total, infected) {
    model$initial_population = list(total = total, infected = infected)
    model
}

##' @family flexmodels
##' @export
add_state_mappings = function(
    model,
    eigen_drop_pattern,
    infected_drop_pattern,
    initial_susceptible_pattern) {

    # TODO: pull state_patterns out of disease_free, because it
    # makes more sense as a general concept -- especially with
    # initial_susceptible_pattern -- and maybe in the future
    # it would be good to put things like vaccination categories
    # in there (OTOH vaccination categories are 'user-defined'
    # concepts whereas initial_susceptible_pattern and
    # infected_drop_pattern are 'general' concepts)
    #model$disease_free$state$patterns =
    model$initialization_mappings =
        list(eigen = eigen_drop_pattern,
             infected = infected_drop_pattern,
             susceptible = initial_susceptible_pattern)
    return(model)
}

# compute indices and pass them to the tmb/c++ side ---------------------

##' Add TMB Indices
##'
##' Add, to a compartmental model, indices used to access appropriate values
##' during simulation and calibration using TMB
##'
##' @param model compartmental model
##' @param another compartmental model with indices for TMB
##' @family flexmodels
##' @export
add_tmb_indices <- function(model) {
    model$tmb_indices <- tmb_indices(model)
    return(model)
}

##' @family flexmodels
##' @export
tmb_indices <- function(model) {
    check_spec_ver_archived()


    if (spec_ver_eq("0.1.0")) {
        sp <- c(model$state, model$params, model$sum_vector)
    } else {
        sp <- c(model$state, model$params)
    }

    indices <- list(make_ratemat_indices = ratemat_indices(model$rates, sp))

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
    if (spec_ver_eq("0.1.1")) {
        indices$disease_free = disease_free_indices(model)
        indices$linearized_params = linearized_param_indices(model)
        indices$outflow = outflow_indices(model$outflow, model$ratemat)
        indices$linearized_outflow = outflow_indices(
            model$linearized_outflow, model$ratemat)
        indices$initialization_mapping = initialization_mapping_indices(model)
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

    ## unpack model structure
    unpack(model)
    unpack(tmb_indices)
    unpack(make_ratemat_indices)
    if (spec_ver_gt("0.0.3")) {
        unpack(timevar$piece_wise)
    }


    ## make ad functions
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
            parameters = list(params = c(params)),
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
            parameters = list(params = c(params)),
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
            parameters = list(params = c(params)),
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
            parameters = list(params = c(params)),
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
            parameters = list(params = c(params)),
            DLL = DLL
        )
    } else if (spec_ver_gt("0.0.6")) { # 0.1.0 and up
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
                linearized_outflow_row_count = linearized_outflow_row_count,
                linearized_outflow_col_count = linearized_outflow_col_count,
                linearized_outflow_rows = linearized_outflow_rows,
                linearized_outflow_cols = linearized_outflow_cols,
                outflow_row_count = outflow_row_count,
                outflow_col_count = outflow_col_count,
                outflow_rows = outflow_rows,
                outflow_cols = outflow_cols,
                do_hazard = do_hazard,
                numIterations = iters
            ),
            parameters = list(params = c(params)),
            DLL = DLL
        )
    } else {
        stop("This feature is not supported by your installation of MacPan")
    }
    return(dd)
}

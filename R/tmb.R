##' Initialize Compartmental Model
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
##' @return object representing a compartmental model
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
                    %>% mutate(tv_spi = find_vec_indices(Symbol, c(state, params)))
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
    }
    if (spec_ver_gt("0.0.4")) model$do_hazard <- do_hazard
    if (spec_ver_gt("0.0.6")) {
        model$sums = list()
        model$sum_vector = numeric(0L)
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

# ----------------------------
# Utilities for rate functions

## regex pattern for finding variables
## (e.g. any parameter or state variable)
## variable_regex looks like this '(beta0|Ca|...|zeta|S|E|Ia|...|V)'
variable_regex <- function(...) {
    return(getOption("MP_name_search_regex"))
    # only works because there are no reserved words allowed yet
    return('[A-z]([A-z][0-9])+')

    # this strategy failed (e.g. Isum matches Is)
    character_class <-
        (list(...)
            %>% lapply(names)
            %>% unlist()
            %>% paste0(collapse = "|")
        )
    paste0("(", character_class, ")", sep = "")
}
get_variables <- function(x) {
    r = gregexpr(variable_regex(), x)
    unlist(regmatches(x, r))
}
## FIXME: this only works because complements (1 - x) and
## inverses (1 / x) are so similar in structure
find_operators <- function(x, operator) {
    grepl(
        paste0("\\( *1 *", operator, " *",
               #variable_regex(params, state),
               variable_regex(),
               collapse = ""
        ), x
    )
}
factor_table <- function(x) {
    data.frame(
        var = unlist(lapply(x, get_variables)),
        compl = unlist(lapply(x, find_operators, "-")),
        invrs = unlist(lapply(x, find_operators, "/"))
    )
}

find_pos_grep <- function(tags, x) {
    pos_list = (tags
     %>% sprintf(fmt = "^%s(_|$)")
     %>% lapply(grep, x)
     %>% unlist(use.names = FALSE)
    )
}

check_in_rate_fun = function(pf) {
    stopifnot(inherits(get("model", envir = pf), 'flexmodel'))
    f = get("formula", envir = pf)
    stopifnot(
        inherits(f, "formula") |
        inherits(f, "struc") |
        inherits(f, "character")
    )
}

#' @export
cross = function(from, to, mat) {
    expand.grid(
        from_pos = unique(find_pos_grep(from, rownames(mat))),
        to_pos = unique(find_pos_grep(to, colnames(mat)))
    )
}

#' @export
pwise = function(from, to, mat) {
    from_pos = find_pos_grep(from, rownames(mat))
    to_pos = find_pos_grep(to, colnames(mat))
    if(length(from_pos) != length(to_pos)) {
        stop("\nargument 'from' matches to ", length(from_pos), " row indices.",
             "\nargument 'to' matches to ", length(to_pos), " column indices.",
             "\nbut these numbers must match for valid pairwise indexing of ",
             "the rate matrix")
    }
    cbind(from_pos, to_pos)
}

#' Paste with Underscore Separator
#' @export
`%_%` = function(x, y) paste(x, y, sep = "_")

#' Paste with Blank Separator
#'
#' Like Python string `+`
#' @export
`%+%` = function(x, y) paste(x, y, sep = "")

# ---------------------
# rate and associated functions:
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

#' @export
mat_rate = function() {
    stop("\nrate specification with matrices is ",
         "coming sometime in the future ... maybe\n",
         "in the meantime you can specify vector-valued rates with vec_rate\n",
         "see this document for potentially more information on priorities:\n",
         options("MP_flex_spec_doc_site")[[1]])
}


#' @export
lookup_pairwise = function(from, to, M) {
    i = pwise(from, to, M)
    data.frame(
        from = rownames(M)[i[,"from_pos"]],
        to = colnames(M)[i[,"to_pos"]]
    )
}


#' Rate Matrix Loopup Table
#'
#' @param state state_pansim object
#' @param ratemat rate matrix
#' @export
rate_matrix_lookup = function(ratemat) {
    ratemat = as(ratemat, "dgTMatrix")
    (data.frame(
        from_pos = ratemat@i + 1,
        to_pos = ratemat@j + 1)
     %>% mutate(
         from_state = ratemat@Dimnames[[1]][from_pos],
         to_state = ratemat@Dimnames[[2]][to_pos]
     )
    )
}

#' Parse a Flexmodel Formula
#'
#' @param x one-sided formula, character vector, or struc object describing
#' the dependence of rate matrix elements on parameter and/or state variables
#' @return depends on the spec version (TODO: add detail once we converge)
#' @export
parse_formula = function(x) {
    y = as.character(x)
    if(inherits(x, "formula")) y = y[[2L]]
    pf = function(y) {
        (y
         %>% strsplit(split = "\\+") %>% getElement(1L)
         %>% strsplit(split = "\\*")
        )
    }
    return(pf(y))
    o = lapply(y, pf)
    if(spec_ver_lt('0.1.0')) {
        return(o[[1L]])
    }
    return(o)
}

find_vec_indices <- function(x, vec) {
    (
        x
            %>% as.character()
            %>% outer(names(vec), "==")
            %>% apply(1, which)
    )
}

#' @export
get_rate_info = function(model, what) lapply(model$rates, '[[', what)

#' @export
get_rate_from = function(model) get_rate_info(model, 'from')

#' @export
get_rate_to = function(model) get_rate_info(model, 'to')

#' @export
get_rate_ratemat_indices = function(model) get_rate_info(model, 'ratemat_indices')

#' @export
get_rate_formula = function(model) get_rate_info(model, 'formula')

#' @export
get_rate_factors = function(model) get_rate_info(model, 'factors')

#' @export
get_rate_state_dependent = function(model) get_rate_info(model, 'state_dependent')

#' @export
get_rate_time_varying = function(model) get_rate_info(model, 'time_varying')

#' @export
get_rate_sum_dependent = function(model) get_rate_info(model, 'sum_dependent')

#' @export
get_n_products = function(model) {
    (model
     %>% get_rate_info('factors')
     %>% lapply('[[', 'prod_indx')
     %>% lapply(length)
    )
}

#' @export
get_n_variables = function(model) {
    (model
     %>% get_rate_info('factors')
     %>% lapply('[[', 'var')
     %>% lapply(length)
    )
}

#' @export
get_n_factors = function(model) {
    (model
     %>% get_rate_info('factors')
     %>% lapply(nrow)
    )
}


#' @export
get_rates = function(model) {
    from = flex_from(model)
    to = flex_to(model)
    mapply(function(from, to) {
        model$ratemat[from, to]
    }, from = from, to = to,
    SIMPLIFY = TRUE)
}

compute_rates = function(model) {
    (model
     %>% get_rate_info("formula")
     %>% lapply(function(x) ifelse(inherits(x, 'formula'), as.character(x[2]), x))
     %>% sapply(as.character)
     %>% struc
     %>% struc_eval(c(as.list(model$params), as.list(model$state), as.list(model$sum_vector)))
     %>% c()
     %>% setNames(names(model$rates))
    )
}

compute_rate_from_indices = function(model, i) {
    unpack(get_indices_per_rate(model, i))
    from = model$tmb_indices$make_ratemat_indices$from[i]
    to = model$tmb_indices$make_ratemat_indices$to[i]
    count = model$tmb_indices$make_ratemat_indices$count[i]
    sp = c(model$state, model$param, model$sum_vector)
    result = 0
    prod = 1.0
    for (j in seq_along(spi)) {
        x = sp[spi[j]]
        if (modifier[j] > 3) {
            result = result + prod
            prod = 1
        }
        if (modifier[j] %in% c(1, 3, 5, 7)) {
            x = 1-x
        } else if (modifier[j] %in% c(2, 3, 6, 7)){
            if (x != 0) {
                x = 1/x
            }
        }
        prod = prod * x
    }
    result = result + prod
    return(result)
}

#' @export
check_rates = function(model, eps = 1e-5) {
    (data.frame(get_rates(test_model), compute_rates(test_model))
        %>% setNames(c('get', 'compute'))
        %>% mutate(diff = abs(get - compute))
        %>% mutate(bads = (diff > eps) | is.nan(compute) | is.na(get))
        %>% filter(bads)
    )
}

#' @param sum name of sum of state variables and parameters
#' @param summands character vector of regular expressions for identifying
#' state variables and parameters to sum together
#' @param state pansim_state object
#' @param params pansim_param object
#' @return indices of summands
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
#' @export
add_state_param_sum = function(model, sum, summands) {
    model$sums[[sum]] = state_param_sum(
        sum, summands, model$state, model$params)

    # assumes that order of sums doesn't change!
    model$sum_vector = get_sum_vector(model$sums)
    model
}

get_sum_vector = function(sums) {
    unlist(lapply(sums, `[[`, 'initial_value'))
}

##' @param x parameter vector
##' @export
has_time_varying <- function(x) {
    spec_check(
        feature = "Time-varying parameters",
        introduced_version = "0.0.3"
    )
    "tv_param_indices" %in% names(attributes(x))
}

##' @export
time_varying_rates <- function(model) {
    model$rates[which_time_varying_rates(model)]
}

##' @export
which_time_varying_rates <- function(model) {
    sd  <- get_rate_info(model, "state_dependent") %>% unlist
    tv  <- get_rate_info(model, "time_varying") %>% unlist
    smd <- get_rate_info(model, "sum_dependent") %>% unlist
    which(sd | tv | smd)
}

##' @export
state_dependent_rates <- function(model) {
    i = get_rate_info(model, "state_dependent") %>% unlist
    model$rates[i]
}

##' @export
sum_dependent_rates = function(model) {
    i = get_rate_info(model, "sum_dependent") %>% unlist
    model$rates[i]
}

##' @export
rate_summary = function(model) {
    data.frame(
        from = get_rate_info(model, "from") %>% unlist,
        to = get_rate_info(model, "to") %>% unlist,
        n_factors = get_n_factors(model) %>% unlist,
        n_products = get_n_products(model) %>% unlist,
        n_variables = get_n_variables(model) %>% unlist,
        state_dependent = get_rate_info(model, "state_dependent") %>% unlist,
        time_varying = get_rate_info(model, "time_varying") %>% unlist,
        sum_dependent = get_rate_info(model, "sum_dependent") %>% unlist
    )
}

##' Add Parallel Accumulators
##'
##' Add parallel accumulators to a compartmental model.
##'
##' @param model TODO
##' @param state_patterns regular expressions for identifying states as
##' parallel accumulators
##' @return another compartmental model with parallel accumulators specified
##' @export
add_parallel_accumulators <- function(model, state_patterns) {
    model$parallel_accumulators <- parallel_accumulators(model, state_patterns)
    return(model)
}

##' @export
parallel_accumulators <- function(model, state_patterns) {
    spec_check(introduced_version = "0.0.2", feature = "Parallel accumulators")
    (
        state_patterns
            %>% lapply(function(x) {
                grep(x, colnames(model$ratemat), value = TRUE)
            })
            %>% unlist()
    )
}

##' Add TMB Indices
##'
##' Add, to a compartmental model, indices used to access appropriate values
##' during simulation and calibration using TMB
##'
##' @param model compartmental model
##' @param another compartmental model with indices for TMB
##' @export
add_tmb_indices <- function(model) {
    model$tmb_indices <- tmb_indices(model)
    return(model)
}

##' @export
sum_indices = function(sums, state, params) {
    sum_index_list = lapply(sums, "[[", "sum_indices")
    sumidx = c(seq_len(length(sums))) + length(state) + length(params)
    sumcount = sapply(sum_index_list, length, USE.NAMES = FALSE) %>% unname
    summandidx = c(unlist(sum_index_list, use.names = FALSE))
    clean_if_empty = function(x) {
        if((length(x) == 0L) | is.null(x)) return(integer(0L))
        x
    }
    list(sumidx = clean_if_empty(sumidx),
         sumcount = clean_if_empty(sumcount),
         summandidx = clean_if_empty(summandidx))
}

##' @export
ratemat_indices <- function(rates, state_params) {
    sp <- state_params
    ratemat_indices <- sapply(rates, `[[`, "ratemat_indices")
    spi <- {
        lapply(rates, function(y) {
            y$factors$var_indx
        }) %>% unlist()
    }
    count <- sapply(rates, function(y) {
        nrow(y$factors)
    })
    modifier <- (rates
        %>% unname
        %>% lapply("[[", "factors")
        %>% bind_rows(.id = "rate_indx")
        %>% mutate(new_rate = as.logical(c(0, diff(as.numeric(rate_indx)))))
        %>% mutate(new_prod = as.logical(c(0, diff(as.numeric(prod_indx)))))
        %>% mutate(add = new_prod & (!new_rate))
        %>% mutate(modifier = 4 * add + 2 * invrs + compl)
        %>% `$`("modifier")
    )
    names(spi) <- colnames(ratemat_indices) <- names(count) <- NULL
    indices <- list(
        from = ratemat_indices[1, ],
        to = ratemat_indices[2, ],
        count = count,
        spi = spi,
        modifier = modifier
    )
    return(indices)
}

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
    return(indices)
}

get_indices_per_rate = function(model, i) {
    unpack(model$tmb_indices$make_ratemat_indices)
    start = sum(c(1, count)[1:i])
    end = start + count[i] - 1L
    list(spi = spi[start:end],
        modifier = modifier[start:end],
        start = start, end = end)
}

##' Make Objective Function with TMB
##'
##' Construct an objective function in TMB from a \code{flexmodel}
##' object. The behaviour of \code{tmb_fun} depends on \code{spec_version()}
##'
##' @param model object of class \code{flexmodel}
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
        if(is.null(timevar$piece_wise)) {
            # HACK: avoid breakage on the C++ side
            #  - proper fix should allow NULLs
            breaks = integer()
            count_of_tv_at_breaks = integer()
            schedule = list()
            schedule$tv_spi = integer()
            schedule$tv_val = numeric()
            schedule$Value = numeric()
            schedule$Type = character()
        } else {
            unpack(timevar$piece_wise)
        }
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
            parameters = list(params = c(params)),
            DLL = DLL
        )
    } else {
        stop("This feature is not supported by your installation of MacPan")
    }
    return(dd)
}


##' Represent a Standard Model as a flexmodel
##'
##' @export
make_unflexmodel <- function(params,
                             state = NULL,
                             start_date = "2020-03-20",
                             end_date = "2020-05-1",
                             params_timevar = NULL,
                             step_args = list()) {

    ## check_start = Sys.time()
    spec_check("0.0.5", "run_sim with TMB")

    ## may need to modify this when we start updating time-varying parameters
    ## on the c++ side (currently params0 should always equal params)
    params0 <- params
    state0 <- state

    step_args_to_flex <- c("do_hazard")
    flex_args <- c(
        list(
            params = params,
            state = state,
            start_date = start_date,
            end_date = end_date,
            params_timevar = params_timevar
        ),
        step_args[names(step_args) %in% step_args_to_flex]
    )
    init_end <- Sys.time()

    ## make_model_start = Sys.time()
    model <- (init_model
        %>% do.call(flex_args)
        %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
        %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
        %>% add_rate("Ia", "R", ~ (gamma_a))
        %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
        %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
        %>% add_rate("Im", "R", ~ (gamma_m))
        %>% add_rate("Is", "H", ~
        (1 - nonhosp_mort) * (phi1) * (gamma_s))
        %>% add_rate("Is", "ICUs", ~
        (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
        %>% add_rate("Is", "ICUd", ~
        (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
        %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
        %>% add_rate("ICUs", "H2", ~ (psi1))
        %>% add_rate("ICUd", "D", ~ (psi2))
        %>% add_rate("H2", "R", ~ (psi3))
        %>% add_rate("H", "R", ~ (rho))
        %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
        %>% add_rate("S", "E", ~
        (Ia) * (beta0) * (1 / N) * (Ca) +
            (Ip) * (beta0) * (1 / N) * (Cp) +
            (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
            (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
        %>% add_parallel_accumulators(c("X", "N", "P", "V"))
        %>% add_tmb_indices()
    )
    return(model)
}

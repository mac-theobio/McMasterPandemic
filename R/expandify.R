## GENERAL UTILITIES
## to expand states with subcategories (age, vaccine status)

## x_y, with x varying faster
expand_names <- function(x, y, sep = "_") {
    unlist(lapply(y, function(a) paste(x, a, sep = sep)))
}

##' distribute counts given a desired distribution (with smart rounding)
##'
##' @param total total count to distribute
##' @param dist target distribution (a vector that sums to 1)
##'
##' @export
distribute_counts <- function(total, dist) {
    ## scale up distribution to total
    scaled_dist <- total * dist
    ## use smart_round in utils to round to whole numbers
    ## while preserving sum of vector
    counts <- smart_round(scaled_dist)
    ## FIXME: smart round always shifts counts to the end of the distribution
    ## (older age groups here)... is this a problem
    return(counts)
}

##' collapse (non-accumulator) states into subcategories (ages, vax status)
##' @param x a state vector or data frame (each row is a different time point)
##' @param values_only just return values (unlisted and unnamed?)
##' @importFrom dplyr matches
##' @export
##' @return a tibble of counts aggregated across epidemiological states
## TODO: rewrite this as a generic function with custom methods for state_pansim
## and pansim objects
condense_state <- function(x, return_type = c("tibble", "named_vector", "unnamed_vector")) {
    return_type <- match.arg(return_type)

    ## R CMD CHECK doesn't understand dplyr, this is a workaround
    obs_number <- subcat <- value <- NULL

    ## define non-accumulator state names
    states <- c(
        "S", "E",
        "I", "H", "ICU",
        "R", "D"
    )
    states_regex <- paste0("^", states)

    ## count number of sub categories
    ## n_subcats <- str_count(names(x[1]), "_")
    ## if(n_subcats == 0) stop("this state vector has not been expanded")
    ## subcat_labels <- paste0("subcat", 1:n_subcats)

    ## if x isn't already a df, convert it to one so we can use tidyverse tools
    if (is.null(nrow(x))) {
        x <- as.data.frame(t(unclass(x)))
    }

    (x
    ## keep only non-accumulator states
        %>% select(matches(paste(states_regex, collapse = "|"),
            ignore.case = FALSE
        ))
        ## append row number to keep value from the same timestep together
        ## before pivoting longer
        ## (sim observations have a date col, but state vectors do not
        ## and i want this function to work in both cases)
        %>% mutate(obs_number = 1:nrow(x))
        ## pivot longer to make state aggregation easier
        %>% pivot_longer(-obs_number,
            names_to = "var"
        )
        ## separate state from subcategories
        %>% separate(var,
            into = c("state", "subcat"),
            sep = "_",
            extra = "merge"
        )
        ## turn subcat col into factor to preserve original ordering
        %>% mutate(subcat = as_factor(subcat))
        ## group by everything except state
        %>% group_by(obs_number, subcat)
        %>% summarise(value = sum(value), .groups = "drop")
        ## put data back into wide format
        %>% pivot_wider(
            names_from = "subcat",
            values_from = "value"
        )
        ## drop observation number
        %>% select(-obs_number)
    ) -> x

    ## repair age cats in names
    x <- repair_names_age(x)

    x <- switch(return_type,
        tibble = x,
        named_vector = tibble_row_to_named_vec(x),
        unnamed_vector = unname(unlist(x))
    )

    return(x)
}

##
## VAXIFY ##
##

##' Construct vector of vaccine category labels
##'
##' @param model_type choose either the one-dose or the two-dose model
##'
##' @export
mk_vaxcats <- function(model_type = c("onedose", "twodose")) {
    model_type <- match.arg(model_type)

    if (model_type == "onedose") {
        cats <- c("unvax", "vaxdose1", "vaxprotect1")
    }

    if (model_type == "twodose") {
        cats <- c("unvax", "vaxdose1", "vaxprotect1", "vaxdose2", "vaxprotect2")
    }

    attr(cats, "model_type") <- model_type

    return(cats)
}

## PARAM TOOLS

##' Edit parameter description attribute to include vax-specific definitions
##' (helper function for `expand_params_vax()`)
##'
##' @param params_desc parameter descriptions, as initialized in `read_params()`
##' @inheritParams mk_vaxcats
##'
##' @export
expand_params_desc_vax <- function(params_desc, model_type = c("onedose", "twodose")) {
    model_type <- match.arg(model_type)

    params_desc[["vax_doses_per_day"]] <- "Total number of doses administered per day in the entire population"
    params_desc[["vax_efficacy_dose1"]] <- "Infection-blocking one-dose efficacy of the vaccine"
    params_desc[["vax_response_rate"]] <- "Average number of days to vaccine-derived immunity, after a dose, for an individual who has never been infected"
    params_desc[["vax_response_rate_R"]] <- "Average number of days to vaccine-derived immunity, after a dose, for an individual who was infected after being dosed (but not yet protected)"
    params_desc[["vax_alpha_dose1"]] <- "Proportion of infections in individuals protected with one dose of the vaccine that are asymptomatic (vs symptomatic)"
    params_desc[["vax_mu_dose1"]] <- "Proportion of infections in individuals protected with one dose of the vaccine that are asymptomatic (vs symptomatic)"

    if (model_type == "twodose") {
        params_desc[["vax_efficacy_dose2"]] <- "Infection-blocking two-dose efficacy of the vaccine"
        params_desc[["vax_alpha_dose2"]] <- "Proportion of infections in individuals protected with two doses of the vaccine that are asymptomatic (vs symptomatic)"
        params_desc[["vax_mu_dose2"]] <- "Proportion of infections in individuals protected with two doses of the vaccine that are asymptomatic (vs symptomatic)"

        params_desc[["vax_prop_first_dose"]] <- "Proportion of the total number of doses administered per day that are first doses"
    }

    return(params_desc)
}

##' Expand parameter list to include vaccination
##'
##' @param params parameter list (e.g. read in with `read_params()`)
##' @inheritParams mk_vaxcats
##' @param vax_doses_per_day total number of vaccine doses administered per day
##' @param vax_prop_first_dose the proportion of vaccine doses administered per day that are first doses
##' @param vax_efficacy_dose1 infection-blocking vaccine efficacy for dose 1
##' @param vax_efficacy_dose2 infection-blocking vaccine efficacy for dose 2
##' @param vax_avg_response_time average number of days it takes for a vaccine dose to confer protection for an individual that has never been infected
##' @param vax_avg_response_time_R average number of days it takes for a vaccine dose to confer protection for an individual that has previously been infected (e.g. between being dosed and eliciting the protective immune response)
##' @param vax_alpha_dose1 proportion of infections that are asymptomatic in individuals that have received one dose of the vaccine
##' @param vax_alpha_dose2 proportion of infections that are asymptomatic in individuals that have received two doses of the vaccine
##' @param vax_mu_dose1 proportion of infections that are mild (vs. severe) in individuals that have received one dose of the vaccine
##' @param vax_mu_dose2 proportion of infections that are mild (vs. severe) in individuals that have received two doses of the vaccine
##'
##' @examples
##' params <- read_params("PHAC.csv")
##' params_vax <- expand_params_vax(params)
##' @export
expand_params_vax <- function(params,
                              model_type = c("onedose", "twodose"),
                              vax_doses_per_day = 1e5,
                              vax_prop_first_dose = 1,
                              vax_efficacy_dose1 = 0.6,
                              vax_efficacy_dose2 = 0.9,
                              vax_avg_response_time = 14,
                              vax_avg_response_time_R = 7,
                              vax_alpha_dose1 = 0.5,
                              vax_alpha_dose2 = 0.5,
                              vax_mu_dose1 = 1,
                              vax_mu_dose2 = 1) {

    ## prep inputs
    ##
    model_type <- match.arg(model_type)
    vax_cat <- mk_vaxcats(model_type)

    ## grab existing parameter descriptions
    params_desc <- attr(params, "description")

    ## convert to list
    params <- as.list(params)

    ## perform updates
    ##
    dose1_pars <- c(
        "vax_doses_per_day",
        "vax_efficacy_dose1",
        "vax_alpha_dose1",
        "vax_mu_dose1"
    )

    for (par in dose1_pars) {
        params[[par]] <- get(par)
    }
    ## update average immune response rate
    params[["vax_response_rate"]] <- 1 / vax_avg_response_time
    params[["vax_response_rate_R"]] <- 1 / vax_avg_response_time_R

    if (model_type == "twodose") {
        dose2_pars <- c(
            "vax_prop_first_dose",
            "vax_efficacy_dose2",
            "vax_alpha_dose2",
            "vax_mu_dose2"
        )

        for (par in dose2_pars) {
            params[[par]] <- get(par)
        }
    }

    ## prep output
    ##

    ## add attributes
    attr(params, "description") <- expand_params_desc_vax(
        params_desc,
        model_type
    )
    attr(params, "vax_cat") <- vax_cat

    ## add class
    class(params) <- "params_pansim"

    return(params)
}

##' condense vaxified parameters to base parameters
##'
##' @param params vaxified parameters (e.g. generated with `expand_params_vax()`)
##'
##' @return an object of class `params_pansim`
##' @export
##'
##' @examples condense_params_vax(expand_params_vax(read_params("PHAC.csv")))
condense_params_vax <- function(params) {
    ## get description of vax params
    desc <- attr(params, "description")

    ## remove vax params and their descriptions
    params <- unlist(params[!grepl("^vax", names(params))])
    desc <- desc[names(params)]

    ## add updated description
    attr(params, "description") <- desc

    ## remove vax cat attribute
    attr(params, "vax_cat") <- NULL

    ## reinstate class
    class(params) <- "params_pansim"

    return(params)
}

## STATE TOOLS

##' expand state vector by vaccination status
##'
##' by default, everyone starts with an unvaccinated status
##'
##' @param x state vector
##' @inheritParams mk_vaxcats
##' @param unif should individuals in each epidemic category be distributed evenly among vaccination strata? (if FALSE, put everyone in the first vaccination category, i.e., the unvaccinated category)
##' @examples
##' params <- read_params("PHAC.csv")
##' ss <- make_state(params=params)
##' ss2 <- expand_state_vax(ss)
##' @export
expand_state_vax <- function(x,
                             model_type = c("onedose", "twodose"),
                             unif = FALSE) {
    model_type <- match.arg(model_type)
    vax_cat <- mk_vaxcats(model_type = model_type)

    ## save attributes
    original_attributes <- attributes(x)

    ## make new names
    new_names <- expand_names(names(x), vax_cat)

    ## set up output
    out <- rep(0, length(new_names))
    names(out) <- new_names

    ## check if original state vector has been normalized to sum to 1
    normalized <- sum(x) == 1

    if (unif) {
        ## distribute people within each epidemic category uniformly across vax strata
        state_regex_list <- paste0("^", attr(x, "epi_cat"))
        ## append age if it exists
        if (has_age(x)) {
            state_regex_list <- expand_names(state_regex_list, sub("\\+", "\\\\+", attr(x, "age_cat")))
        }

        for (state_regex in state_regex_list) {
            if (normalized) {
                out[grepl(paste0(state_regex, "_"), names(out))] <- rep(x[grepl(paste0(state_regex, "$"), names(x))], length(vax_cat)) / length(vax_cat)
            } else {
                out[grepl(paste0(state_regex, "_"), names(out))] <- distribute_counts(
                    total = x[grepl(paste0(state_regex, "$"), names(x))],
                    dist = rep(1, length(vax_cat)) / length(vax_cat)
                )
            }
        }
    } else {
        ## put everyone in the unvaccinated class by default (first category)
        out[grepl(vax_cat[1], names(out))] <- x
    }

    ## update names in original attributes to new names
    original_attributes$names <- new_names

    ## restore attributes to output
    attributes(out) <- original_attributes

    ## add vax cat attribute to output
    attr(out, "vax_cat") <- vax_cat

    return(out)
}

## CONDENSE STATE OR SIM RESULTS

##' collapse vaccination categories (only; don't do other condensation)
##' @param x a state vector or simulation result df (an object of class `state_pansim` or `pansim`)
##'
##' @importFrom stringr str_count
##' @importFrom tidyr unite
##' @importFrom dplyr all_of
##' @importFrom dplyr everything
##' @export
## wrap into condense() ?
condense_vax <- function(x) {
    ## R CMD CHECK doesn't understand dplyr, this is a workaround
    state <- value <- subcat1 <- . <- NULL
    ## save original attributes
    original_attributes <- attributes(x)

    ## x has age?
    age <- has_age(x)

    ## get input type
    input_class <- class(x)

    ## figure out how many types of subcategories there are
    n_subcats <- str_count(names(x)[2], "_")
    subcat_labels <- paste0("subcat", 1:n_subcats)

    ## if x is a vector (e.g. a state vec), convert it to a (wide) df
    if (is.null(nrow(x))) {
        x <- as.data.frame(t(unclass(x)))
    }

    ## check if there are columns for time
    time_vars <- length(intersect(c("t", "date"), names(x))) > 0

    ## pivot longer to make state aggregation easier
    ## take care to preserve t and/or date columns, if they exist
    if (time_vars) {
        time_names <- intersect(c("t", "date"), names(x))
        (x
            %>% pivot_longer(-all_of(time_names),
                names_to = "var"
            )
        ) -> x
    } else {
        (x
            %>% pivot_longer(everything(),
                names_to = "var"
            )
        ) -> x
    }

    (x
    ## separate state from subcategory
        %>% separate(var,
            into = c("state", subcat_labels),
            sep = "_",
            extra = "merge"
        )
        ## turn state column into factor to keep original ordering
        %>% mutate(state = factor(state, levels = unique(state)))
    ) -> x

    ## FIXME: figure out which subcategory includes vaccination
    ## and then group by all but that column (also except values, of course)
    ## currently, this is hacky

    ## aggregate value by state (and timestep, if it exists)
    group_names <- intersect(
        c("t", "date", "state"),
        names(x)
    )
    if (age) group_names <- c(group_names, "subcat1")

    (x
        %>% group_by(across(group_names))
    ) -> x

    (x
        %>% summarise(value = sum(value), .groups = "drop")
        %>%
        {
            if (age) {
                unite(
                    ., "state",
                    state, subcat1
                )
            } else {
                .
            }
        }
        ## put data back into wide format
        %>% pivot_wider(
            names_from = "state",
            values_from = "value"
        )
    ) -> x

    ## finalize output type
    if ("state_pansim" %in% input_class) {
        x <- unlist(x, use.names = TRUE)
    }

    if ("pansim" %in% input_class) {
        x <- as.data.frame(x)
    }

    ## repair age name (if present); otherwise, do nothing
    x <- repair_names_age(x)

    ## add back original attributes
    new_names <- names(x)
    original_attributes$names <- new_names
    attributes(x) <- original_attributes

    ## remove vax categories attributes
    attr(x, "vax_cat") <- NULL

    return(x)
}

## SIM TOOLS

##' generate per capita daily vaccination rates
##' @param state state vector (an object of class `state_pansim`)
##' @param params model parameters (an object of class `params_pansim`)
##' @export
make_vaxrate <- function(state, params) {

    ## initialize output vec
    vax_rate <- c()

    if (!has_vax(state) | !has_vax(params)) stop("need vaxified state and params to make vaccination rates")

    model_type <- get_vax_model_type(get_vax(params))

    ## pull out non-symptomatic *and* unvaccinated states
    asymp_unvax_regex <- sprintf(
        "^(%s)_.*unvax",
        paste(asymp_cat, collapse = "|")
    )

    ## same as above but for pop that is protected by first dose
    if (model_type == "twodose") {
        asymp_vaxprotect1_regex <- sprintf(
            "^(%s)_.*vaxprotect1",
            paste(asymp_cat, collapse = "|")
        )
    }

    asymp_unvax_N <- condense_state(
        state[grepl(asymp_unvax_regex, names(state))],
        return_type = "named_vector"
    )
    ## same as above but for pop that is protected by first dose
    if (model_type == "twodose") {
        asymp_vaxprotect1_N <- condense_state(
            state[grepl(asymp_vaxprotect1_regex, names(state))],
            return_type = "named_vector"
        )
    }

    if (model_type == "onedose") {
        x <- params[["vax_doses_per_day"]] / asymp_unvax_N
        x[is.nan(x)] <- 0 ## replace NaN with 0, which occurs when asymp_unvax_N is 0
        vax_rate$dose1 <- x
    }

    if (model_type == "twodose") {
        x <- params[["vax_prop_first_dose"]] * params[["vax_doses_per_day"]] / asymp_unvax_N
        x[is.nan(x)] <- 0 ## replace NaN with 0, which occurs when asymp_unvax_N is 0
        vax_rate$dose1 <- x

        x <- (1 - params[["vax_prop_first_dose"]]) * params[["vax_doses_per_day"]] / asymp_vaxprotect1_N
        x[is.nan(x)] <- 0 ## replace NaN with 0, which occurs when asymp_unvax_N is 0
        vax_rate$dose2 <- x
    }

    return(vax_rate)
}

##' compute and add vaxrates to ratemat
##' (returns whole ratemat because the update part is non-trivial)
##' @param state state vector (an object of class `state_pansim`)
##' @param params model parameters (an object of class `params_pansim`)
##' @param ratemat model rate matrix
##' @export
add_updated_vaxrate <- function(state, params, ratemat) {
    vax_cat <- get_vax(params)
    model_type <- get_vax_model_type(vax_cat)

    ## if original matrix isn't sparse, make it sparse for this bit
    original_sparse <- inherits(ratemat, "Matrix")
    if (!original_sparse) {
        ratemat <- Matrix::Matrix(ratemat)
    }

    ## calculate per capita rate of doses per day
    ## (per capita = per non-symptomatic individuals here, because that's
    ## who's getting vaccinated)
    vax_rate <- make_vaxrate(state, params)

    epi_states <- attr(state, "epi_cat")

    ## update unvax -> vaxdose1 block (& vaxprotect1 -> vaxdose2 block)
    if (!has_age(params)) {
        ## set up block diagonal matrix for vaccine allocation step within each age group
        vax_block_dose1 <- matrix(0,
            nrow = length(epi_states),
            ncol = length(epi_states),
            dimnames = list(epi_states, epi_states)
        )
        if (model_type == "twodose") vax_block_dose2 <- vax_block_dose1

        ## for every epi state getting vaccinated (non-symptomatic states), assign vax rate between matching epi states, and add rate for flow into vax accumulator compartment
        for (state_cat in asymp_cat) {
            ## unvax to vaxdose1
            index <- pfun(
                paste0(state_cat),
                paste0(state_cat),
                vax_block_dose1
            )
            vax_block_dose1[index] <- vax_rate$dose1
            if (model_type == "twodose") vax_block_dose2[index] <- vax_rate$dose2

            ## accumulator (not present e.g. when we do rExp in make_state)
            if ("V" %in% epi_states) {
                index <- pfun(
                    paste0(state_cat),
                    paste0("V"),
                    vax_block_dose1
                )
                vax_block_dose1[index] <- vax_rate$dose1
                if (model_type == "twodose") vax_block_dose2[index] <- vax_rate$dose2
            }
        }

        ## convert vax_block to Matrix::Matrix object for subset assignement
        vax_block_dose1 <- Matrix::Matrix(vax_block_dose1)
        if (model_type == "twodose") vax_block_dose2 <- Matrix::Matrix(vax_block_dose2)

        ## just once, without ages
        from_regex <- vax_cat[1] ## unvax
        to_regex <- vax_cat[2]
        ratemat[
            grepl(from_regex, dimnames(ratemat)$from),
            grepl(to_regex, dimnames(ratemat)$to)
        ] <- vax_block_dose1
        if (model_type == "twodose") {
            ## just once, without ages
            from_regex <- vax_cat[3] ## vaxprotect1
            to_regex <- vax_cat[4] ## vaxdose2
            ratemat[
                grepl(from_regex, dimnames(ratemat)$from),
                grepl(to_regex, dimnames(ratemat)$to)
            ] <- vax_block_dose2
        }
    } else {
        ## for each age
        for (age in attr(params, "age_cat")) {
            from_regex_suffix <- sub(
                "\\+", "\\\\+",
                paste0(age, "_", vax_cat[1])
            )
            to_regex_suffix <- sub(
                "\\+", "\\\\+",
                paste0(age, "_", vax_cat[2])
            )

            dose1_rate <- vax_rate$dose1[grepl(from_regex_suffix, names(vax_rate$dose1))]

            from_regex <- paste0("^(S|E|Ia|Ip|R)_", from_regex_suffix)
            to_regex <- paste0("^(S|E|Ia|Ip|R)_", to_regex_suffix)
            block_size <- 5

            ratemat[
                grepl(from_regex, dimnames(ratemat)$from),
                grepl(to_regex, dimnames(ratemat)$to)
            ] <- diag(dose1_rate, nrow = block_size, ncol = block_size, names = FALSE)

            ## add flows into vax accumulator compartment
            if ("V" %in% epi_states) {
                ratemat[
                    grepl(from_regex, dimnames(ratemat)$from),
                    grepl(paste0("^V_", to_regex_suffix), dimnames(ratemat)$to)
                ] <- rep(as.numeric(dose1_rate), block_size)
            }

            if (model_type == "twodose") {
                from_regex_suffix <- sub(
                    "\\+", "\\\\+",
                    paste0(age, "_", vax_cat[3])
                )
                to_regex_suffix <- sub(
                    "\\+", "\\\\+",
                    paste0(age, "_", vax_cat[4])
                )

                dose2_rate <- vax_rate$dose2[grepl(from_regex_suffix, names(vax_rate$dose2))]

                from_regex <- paste0("^(S|E|Ia|Ip|R)_", from_regex_suffix)
                to_regex <- paste0("^(S|E|Ia|Ip|R)_", to_regex_suffix)
                block_size <- 5

                ratemat[
                    grepl(from_regex, dimnames(ratemat)$from),
                    grepl(to_regex, dimnames(ratemat)$to)
                ] <- diag(dose2_rate, nrow = block_size, ncol = block_size, names = FALSE)

                if ("V" %in% epi_states) {
                    ratemat[
                        grepl(from_regex, dimnames(ratemat)$from),
                        grepl(paste0("^V_", to_regex_suffix), dimnames(ratemat)$to)
                    ] <- rep(as.numeric(dose2_rate), block_size)
                }
            }
        }
    }

    ## check that calculated per capita vax rate per day squares with total number of daily doses specified in params
    if (model_type == "onedose") {
        ratemat_dose1_subset <- ratemat[
            grepl(vax_cat[1], dimnames(ratemat)$from),
            grepl(vax_cat[2], dimnames(ratemat)$to)
        ]
        state_dose1_subset <- state[grepl(vax_cat[1], names(state))]
        ratemat_doses <- sum(ratemat_dose1_subset %*% state_dose1_subset)
        total_doses_match <- isTRUE(all.equal(ratemat_doses, params[["vax_doses_per_day"]]))
        ## if total doses allocated via rate matrix does not match match total doses specified in params
        if (!total_doses_match) {
            ## *and* it's not because we've depleted the population eligible for vaccination
            if (sum(state_dose1_subset) >= params[["vax_doses_per_day"]]) warning("calculated daily vax rate exceeds size of remaining population eligible for vaccination; make sure you're using do_hazard = TRUE to avoid jumps into negative state variables")
        }
    }

    if (model_type == "twodose") {
        ## check dose 1
        ratemat_dose1_subset <- ratemat[
            grepl(vax_cat[1], dimnames(ratemat)$from),
            grepl(vax_cat[2], dimnames(ratemat)$to)
        ]
        state_dose1_subset <- state[grepl(vax_cat[1], names(state))]
        ratemat_dose1 <- sum(ratemat_dose1_subset %*% state_dose1_subset)
        total_dose1_match <- isTRUE(all.equal(ratemat_dose1, sum(params[["vax_prop_first_dose"]] * params[["vax_doses_per_day"]])))
        ## if total doses allocated via rate matrix does not match match total doses specified in params
        if (!total_dose1_match) {
            ## *and* it's not because we've depleted the population eligible for vaccination
            if (sum(state_dose1_subset) >= sum(params[["vax_prop_first_dose"]] * params[["vax_doses_per_day"]])) warning("calculated daily vax dose 1 rate exceeds size of remaining population eligible for dose 1; make sure you're using do_hazard = TRUE to avoid jumps into negative state variables")
        }

        ## check dose 2
        ratemat_dose2_subset <- ratemat[
            grepl(vax_cat[3], dimnames(ratemat)$from),
            grepl(vax_cat[4], dimnames(ratemat)$to)
        ]
        state_dose2_subset <- state[grepl(vax_cat[3], names(state))]
        ratemat_dose2 <- sum(ratemat_dose2_subset %*% state_dose2_subset)
        total_dose2_match <- isTRUE(all.equal(ratemat_dose2, sum((1 - params[["vax_prop_first_dose"]]) * params[["vax_doses_per_day"]])))
        ## if total doses allocated via rate matrix does not match match total doses specified in params
        if (!total_dose2_match) {
            ## *and* it's not because we've depleted the population eligible for vaccination
            if (sum(state_dose2_subset) >= sum((1 - params[["vax_prop_first_dose"]]) * params[["vax_doses_per_day"]])) warning("calculated daily vax dose 2 rate exceeds size of remaining population eligible for dose 2; make sure you're using do_hazard = TRUE to avoid jumps into negative state variables")
        }
    }


    ## make updated ratemat have the same type and attributes as original ratemat
    if (!original_sparse) {
        ratemat <- as.matrix(ratemat)
    }

    return(ratemat)
}

##' Compute total number of vaccine doses actually administered per day in a simulation
##'
##' @param res simulation result from model with vaccination (generated using `run_sim`)
##' @param dose the number of doses to consider, vector. (note that Moderna/Pfizer both recommend 2 doses) (??)
##'
##' @return `tibble` with columns for date and total vaccine doses actually administered each day (as captured in the simulation)
##' @export
get_doses_per_day <- function(res, dose = c(1, 2)) {
    dose <- match.arg(dose)

    if (!has_vax(res)) {
        stop("simulation result must be from model with vaccination")
    }

    ## get vax cat that dosed individuals flow into, depending on which dose
    ## for dose 1, its the 2nd category (vaxdose1),
    ## for dose 2, it's the 4th category (vaxdose2)
    vax_cat <- attr(attr(res, "params"), "vax_cat")[ifelse(dose == 1, 2,
        4
    )]

    print(paste0("computing vax rate based on timeseries for V_", vax_cat, " accumulator compartment"))
    difference <- diff(res[, paste0("V_", vax_cat)])
    date <- res$date[1:(nrow(res) - 1)]

    out <- data.frame(
        date = date,
        total_doses_per_day = difference
    )

    return(out)
}

##
## VARIANTIFY ##
##

## PARAM TOOLS

##' Edit parameter description attribute to include variant-specific definitions
##' (helper function for `expand_params_variant()`)
##'
##' @param params_desc parameter descriptions, as initialized in `read_params()`
##'
##' @export
expand_params_desc_variant <- function(params_desc) {
    params_desc[["variant_prop"]] <- "Current proportion of infections caused by the variant"
    params_desc[["variant_advantage"]] <- "Transmissibility advantage of variant compared to current dominant strain (as a multiplicative factor), e.g., a 50% more transmissible variant would have variant_advantage = 1.5."
    params_desc[["variant_vax_efficacy_dose1"]] <- "One-dose vaccine efficacy against variant (only used in vaxified model)"
    params_desc[["variant_vax_efficacy_dose2"]] <- "Two-dose vaccine efficacy against variant (only used in vaxified model)"

    return(params_desc)
}

##' Expand parameter list to include variant strain
##'
##' @param params parameter list (e.g. read in with \code{read_params()})
##' @param variant_prop current proportion of infections caused by the variant
##' @param variant_advantage transmissibility advantage of variant compared to current dominant strain (as a multiplicative factor), e.g., a 50\% more transmissible variant would have \code{variant_advantage = 1.5}.
##' @param variant_vax_efficacy_dose1 one-dose vaccine efficacy against variant (only used in vaxified model)
##' @param variant_vax_efficacy_dose2 two-dose vaccine efficacy against variant (only used in vaxified model)
##'
##' @examples
##' params <- read_params("PHAC.csv")
##' params_variant <- expand_params_variant(params)
##' @export
expand_params_variant <- function(params,
                                  variant_prop = 0.01,
                                  variant_advantage = 1.5,
                                  variant_vax_efficacy_dose1 = 0.3,
                                  variant_vax_efficacy_dose2 = 0.8) {

    ## grab existing parameter descriptions
    params_desc <- attr(params, "description")

    ## convert to list
    params <- as.list(params)

    ## perform updates
    ##
    par_list <- c(
        "variant_prop",
        "variant_advantage",
        "variant_vax_efficacy_dose1",
        "variant_vax_efficacy_dose2"
    )

    for (par in par_list) {
        params[[par]] <- get(par)
    }

    ## prep output
    ##

    ## add attributes
    attr(params, "description") <- expand_params_desc_variant(params_desc)
    attr(params, "do_variant") <- TRUE

    ## add class
    class(params) <- "params_pansim"

    return(params)
}

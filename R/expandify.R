## GENERAL UTILITIES
## to expand states with subcategories (age, vaccine status)

## x_y, with x varying faster
expand_names <- function(x,y,sep="_") {
  unlist(lapply(y, function(a) paste(x, a, sep=sep)))
}

#' distribute counts given a desired distribution (with smart rounding)
#'
#' @param total total count to distribute
#' @param dist target distribution (a vector that sums to 1)
#' @export
distribute_counts <- function(total, dist){
  ## scale up distribution to total
  scaled_dist <- total*dist
  ## use smart_round in utils to round to whole numbers
  ## while preserving sum of vector
  counts <- smart_round(scaled_dist)
  ## FIXME: smart round always shifts counts to the end of the distribution
  ## (older age groups here)... is this a problem
  return(counts)
}

#' collapse (non-accumulator) states into subcategories (ages, vax status)
#' @param x a state vector or data frame (each row is a different time point)
#' @param values_only just return values (unlisted and unnamed?)
#' @export
#' @return a tibble of counts aggregated across epidemiological states
## TODO: rewrite this as a generic function with custom methods for state_pansim
## and pansim objects
condense_state <- function(x, values_only = FALSE){

  ## define non-accumulator state names
  states <- c("S", "E",
              "I", "H", "ICU",
              "R", "D")
  states_regex <- paste0("^", states)

  ## count number of sub categories
  # n_subcats <- str_count(names(x[1]), "_")
  # if(n_subcats == 0) stop("this state vector has not been expanded")
  # subcat_labels <- paste0("subcat", 1:n_subcats)

  ## if x isn't already a df, convert it to one so we can use tidyverse tools
  if(is.null(nrow(x))){
    x <- as.data.frame(t(unclass(x)))
  }

  (x
    ## keep only non-accumulator states
    %>% select(matches(paste(states_regex, collapse = "|"),
                       ignore.case = FALSE))
    ## append row number to keep value from the same timestep together
    ## before pivoting longer
    ## (sim observations have a date col, but state vectors do not
    ## and i want this function to work in both cases)
    %>% mutate(obs_number = 1:nrow(x))
    ## pivot longer to make state aggregation easier
    %>% pivot_longer(-obs_number,
                     names_to = "var")
    ## separate state from subcategories
    %>% separate(var,
                 into = c("state", "subcat"),
                 sep = "_",
                 extra = "merge")
    ## turn subcat col into factor to preserve original ordering
    %>% mutate(subcat = as.factor(subcat))
    ## group by everything except state
    %>% group_by(obs_number, subcat)
    %>% summarise(value = sum(value), .groups = "drop")
    ## put data back into wide format
    %>% pivot_wider(names_from = "subcat",
                    values_from = "value")
    ## drop observation number
    %>% select(-obs_number)
  ) -> x

  ## repair age cats in names
  x <- repair_names_age(x)

  if(values_only) x <- unname(unlist(x))

  return(x)
}

## VAXIFY

#' Construct vector of vaccine category labels
#'
#' @param use_doses (logical) should we make separate categories for two doses?
#' @export
mk_vaxcats <- function(use_doses = FALSE) {
  if (use_doses) return(c("unvax", "onevax", "twovax"))
  return(c("unvax", "vax"))
}

#' expand state vector by vaccination status
#'
#' by default, everyone starts with an unvaccinated status
#'
#' @param x state vector
#' @param vax_cat vaccine status categories
#' @examples
#' params <- read_params("PHAC_testify.csv")
#' ss <- make_state(params=params)
#' ss2 <- expand_state_vax(ss)
#' @export
## FIXME: make it so that we can start a simulation part-way, with some vaccination
expand_state_vax <- function(x, vax_cat = mk_vaxcats()) {

  ## save attributes
  original_attributes <- attributes(x)

  ## make new names
  new_names <- expand_names(names(x), vax_cat)

  ## set up output
  out <- rep(0, length(new_names))
  names(out) <- new_names

  ## put everyone in the unvaccinated class by default (first category)
  out[grepl(vax_cat[1], names(out))] <- x

  ## update names in original attributes to new names
  original_attributes$names <- new_names

  ## restore attributes to output
  attributes(out) <- original_attributes

  ## add vax cat attribute to output
  attr(out, "vax_cat") <- vax_cat

  return(out)
}

#' collapse vaccination categories (only; don't do other condensation)
#' @param x a state vector or simulation result df (an object of class `state_pansim` or `pansim`)
#' @importFrom stringr str_count
#' @export
## wrap into condense() ?
condense_vax <- function(x) {

  ## get input type
  input_class <- class(x)

  ## figure out how many types of subcategories there are
  n_subcats <- str_count(names(x)[1], "_")
  subcat_labels <- paste0("subcat", 1:n_subcats)

  ## if x is a vector (e.g. a state vec), convert it to a (wide) df
  if(is.null(nrow(x))){
    x <- as.data.frame(t(unclass(x)))
  }

  ## pivot longer to make state aggregation easier
  ## take care to preserve date column, if it exists
  if ("pansim" %in% input_class){
    (x
     %>% pivot_longer(-date,
                      names_to = "var")
    ) -> x
  } else {
    (x
     %>% pivot_longer(everything(),
                      names_to = "var")
    ) -> x
  }

  (x
    ## separate state from subcategory
    %>% separate(var,
                 into = c("state", subcat_labels),
                 sep = "_",
                 extra = "merge")
    ## turn state column into factor to keep original ordering
    %>% mutate(state = factor(state, levels = unique(state)))
  ) -> x

  ## FIXME: figure out which subcategory includes vaccination
  ## and then group by all but that column (also except values, of course)

  ## aggregate value by state (and timestep, if it exists)
  if("pansim" %in% input_class){
    (x
     %>% group_by(date, state)
    ) -> x
  } else {
    (x
     %>% group_by(state)
    ) -> x
  }

  (x
    %>% summarise(value = sum(value), .groups = "drop")
    ## put data back into wide format
    %>% pivot_wider(names_from = "state",
                    values_from = "value")
  )  -> x

  # finalize output type
  if("state_pansim" %in% input_class){
    x <- unlist(x, use.names = TRUE)
  }

  if("pansim" %in% input_class){
    x <- as.data.frame(x)
  }

  class(x) <- input_class

  return(x)
}

## FIXME: write a function that gets doses_per_day vector from data
## with names compatible with use in mk_vaxrates()


#' generate per capita daily vaccination rates
#' @param params model parameters (an object of class `params_pansim`)
#' @param doses_per_day named list of doses per day in the region by vax category and any other subcat (e.g. age)
#' @export
mk_vaxrates <- function(params,
                        doses_per_day = list(vax = 10000)
){

  ## make rate names, depending on whether or not we have age-structured params
  if(has_age(params)){
    rate_names <- expand_names(attr(params, "age_cat"), names(doses_per_day))
  } else {
    rate_names <- names(doses_per_day)
  }

  ## set up vax_rate vector
  vax_rate <- rep(0, length(rate_names))
  names(vax_rate) <- rate_names

  ## updated rates for compartments where there is actually vaccination
  ## FIXME: check if this works properly for ageified parameters
  for(this_name in names(doses_per_day)){
    vax_rate[grepl(paste0("^?_?", this_name),
                   names(vax_rate))] <- doses_per_day[[grepl(paste0("^?_?", this_name),
                                                            names(doses_per_day))]]/params[["N"]]
  }

  return(vax_rate)
}

#' Edit parameter description attribute to include vax-specific definitions
#' (helper function for `expand_params_vax()`)
#'
#' @param params_desc parameter descriptions, as initialized in `read_params()`
#'
#' @return
#' @export
expand_params_desc_vax <- function(params_desc){

  params_desc[["vax_rate"]] <- "Per capita daily vaccination rate, by vaccination category"
  params_desc[["avg_time_to_vax_immunity"]] <- "Average number of days to vaccine-derived immunity, after dose"

  return(params_desc)
}

#' Expand parameter list to include vaccination
#'
#' @param params parameter list (e.g. read in with `read_params()`)
#' @param vax_cat vector of vaccination categories
#' @examples
#' params <- read_params("PHAC_testify.csv")
#' params_age <- expand_params_age(params)
#' @export
expand_params_vax <- function(params,
                              vax_cat = mk_vaxcats(),
                              doses_per_day = NULL){
  ## prep inputs
  #################

  ## grab existing parameter descriptions
  params_desc <- attr(params, "description")

  ## convert to list
  params <- as.list(params)

  ## make vaxrates
  vaxrates_args <- list(params = params)

  if(!is.null(doses_per_day)) vaxrates_args[["doses_per_day"]] <- doses_per_day

  vax_rates <- do.call(mk_vaxrates, vaxrates_args)

  ## perform updates
  #####################

  ## update vax rates (S_unvax -> S_vax)
  params[["vax_rates"]] <- vax_rates
  ## update time to immunity
  params[["avg_time_to_vax_immunity"]] <- 1/14

  ## prep output
  ################

  ## add attributes
  attr(params, "description") <- expand_params_desc_vax(params_desc)
  attr(params, "vax_cat") <- vax_cat

  ## add class
  class(params) <- "params_pansim"

  return(params)
}

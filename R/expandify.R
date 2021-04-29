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
#' @export
mk_vaxcats <- function() {
  return(c("unvax", "vaxdose1", "vaxprotect1"))
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
  n_subcats <- str_count(names(x)[2], "_")
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
#' @param state state vector (an object of class `state_pansim`)
#' @param params model parameters (an object of class `params_pansim`)
#' @export
make_vaxrate <- function(state, params){

  if(!has_vax(state) | !has_vax(params)) stop("need vaxified state and params to make vaccination rates")

  ## pull out non-symptomatic *and* unvaccinated states
  asymp_unvax_regex <- sprintf("^(%s)_.*unvax",
                               paste(asymp_cat, collapse="|"))

  ## FIXME: get this working for age-specific vax_doses_per_day
  ## don't sum over all ages, keep pop-size separate for each age
  asymp_unvax_N <- rowSums(condense_state(
    state[grepl(asymp_unvax_regex, names(state))]
    ))

  ## should be a scalar if we're not doing age-specific stuff
  vax_rate <- params[["vax_doses_per_day"]]/asymp_unvax_N

  return(vax_rate)
}

#' compute and add vaxrates to ratemat
#' (returns whole ratemat because the update part is non-trivial)
#' @param state state vector (an object of class `state_pansim`)
#' @param params model parameters (an object of class `params_pansim`)
#' @param ratemat model rate matrix
#' @export
add_updated_vaxrate <- function(state, params, ratemat){
  vax_cat <- get_vax(params)

  ## capture initial state of ratemat
  if (inherits(ratemat,"Matrix")) {
    sparse <- TRUE
    saved_attrs <- setNames(lapply(aa,attr,x=ratemat),aa)
  } else {
    sparse <- FALSE
  }

  ## convert to a Matrix::Matrix object, so we can assign to subsets within the matrix
  ratemat <- Matrix::Matrix(ratemat)

  ## calculate per capita rate of doses per day
  ## (per capita = per non-symptomatic individuals here, because that's
  ## who's getting vaccinated)
  vax_rate <- make_vaxrate(state, params)

  ## set up block diagonal matrix for vaccine allocation step within each age group
  epi_states <- names(condense_age(condense_vax(state)))
  vax_block <- matrix(0,
                      nrow = length(epi_states),
                      ncol = length(epi_states),
                      dimnames = list(epi_states, epi_states))

  ## for every epi state getting vaccinated (non-symptomatic states), assign vax rate between matching epi states,
  for(state_cat in asymp_cat){
    index <- pfun(paste0(state_cat),
                  paste0(state_cat),
                  vax_block)
    vax_block[index] <- vax_rate
  }

  ## convert vax_block to Matrix::Matrix object for subset assignement
  vax_block <- Matrix::Matrix(vax_block)

  ## update unvax -> vaxdose block
  if(!has_age(params)){
    ## just once, without ages
    from_regex <- vax_cat[1]
    to_regex <- vax_cat[2]
    ratemat[grepl(from_regex, dimnames(ratemat)$from),
      grepl(to_regex, dimnames(ratemat)$to)] <- vax_block
  } else {
    ## for each age
    for(age in attr(params, "age_cat")){
      from_regex <- sub("\\+", "\\\\+",
                        paste0(age, "_", vax_cat[1]))
      to_regex <- sub("\\+", "\\\\+",
                      paste0(age, "_", vax_cat[2]))
      ratemat[grepl(from_regex, dimnames(ratemat)$from),
        grepl(to_regex, dimnames(ratemat)$to)] <- vax_block
    }
  }

  ## check that calculated per capita vax rate per day squares with total number of daily doses specified in params
  ratemat_vax_subset <- ratemat[
    grepl(vax_cat[1], dimnames(ratemat)$from),
    grepl(vax_cat[2], dimnames(ratemat)$to)]
  state_vax_subset <- state[grepl(vax_cat[1], names(state))]
  ratemat_doses <- sum(ratemat_vax_subset %*% state_vax_subset)
  total_doses_match <- all.equal(ratemat_doses, params[["vax_doses_per_day"]])
  ## if total doses allocated via rate matrix does not match match total doses specified in params
  if(!total_doses_match){
    ## *and* it's not because we've depleted the population eligible for vaccination
    if(sum(state_vax_subset)>=params[["vax_doses_per_day"]]) stop("calculated daily vax rate exceeds size of remaining population eligible for vaccination")
  }

  ## make updated ratemat have the same type and attributes as original ratemat
  if (sparse){
    for (a in aa) {
      attr(ratemat,a) <- saved_attrs[[a]]
    }
  } else {
    ratemat <- as.matrix(ratemat)
  }

  return(ratemat)
}

#' Edit parameter description attribute to include vax-specific definitions
#' (helper function for `expand_params_vax()`)
#'
#' @param params_desc parameter descriptions, as initialized in `read_params()`
#'
#' @return
#' @export
expand_params_desc_vax <- function(params_desc){

  params_desc[["vax_doses_per_day"]] <- "Total number of doses administered per day in the entire population"
  params_desc[["vax_efficacy"]] <- "Infection-blocking efficacy of the vaccine"
  params_desc[["vax_response_rate"]] <- "Average number of days to vaccine-derived immunity, after dose"

  return(params_desc)
}

#' Expand parameter list to include vaccination
#'
#' @param params parameter list (e.g. read in with `read_params()`)
#' @param vax_cat vector of vaccination categories
#' @param doses_per_date total number of doses administered per day
#' @param vax_efficacy efficacy of first (only) dose (currently)
#' @param vax_alpha alpha (proportion of infections that are asymptomatic) for individuals who are vaccinated and have had their immune response
#' @param vax_mu mu (proportion of symptomatic infections that are mild) for individuals who are vaccinated and have had their immune response
#'
#' @examples
#' params <- read_params("PHAC_testify.csv")
#' params_age <- expand_params_age(params)
#' @export
expand_params_vax <- function(params,
                              vax_cat = mk_vaxcats(),
                              vax_doses_per_day = 1e5,
                              vax_efficacy = 0.7,
                              vax_avg_response_time  = 14,
                              vax_alpha = 0.5,
                              vax_mu = 1){
  ## prep inputs
  #################

  ## grab existing parameter descriptions
  params_desc <- attr(params, "description")

  ## convert to list
  params <- as.list(params)

  ## perform updates
  #####################

  ## add doses per day
  params[["vax_doses_per_day"]] <- vax_doses_per_day
  ## add efficacy
  params[["vax_efficacy"]] <- vax_efficacy
  ## update average immune response rate
  params[["vax_response_rate"]] <- 1/vax_avg_response_time
  ## add updates to epi parameters for vaxprotect layer
  params[["vax_alpha"]] <- vax_alpha
  params[["vax_mu"]] <- vax_mu

  ## prep output
  ################

  ## add attributes
  attr(params, "description") <- expand_params_desc_vax(params_desc)
  attr(params, "vax_cat") <- vax_cat

  ## add class
  class(params) <- "params_pansim"

  return(params)
}

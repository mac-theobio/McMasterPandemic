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


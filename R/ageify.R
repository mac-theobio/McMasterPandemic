#' construct age categories
#' @param min minimum age
#' @param max maximum age
#' @param da age bin width
#' @export
mk_agecats <- function(min=0,max=90,da=10) {
  s1 <- seq(min,max,by=da)
  ## if we have one-year age groups, don't make hyphenated
  ## categories
  if(da == 1){
    out <- c(as.character(s1[-length(s1)]),
             paste0(s1[length(s1)], "+"))
  } else {
    out <- c(sprintf("%d-%d",s1[-length(s1)],(s1[-1]-1)),
             paste0(s1[length(s1)],"+"))
  }
  return(out)
}

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

## STATES

#' expand state vector and rate matrix by age classes and population distribution
#'
#' epidemiological state varies fast, age category varies slowly
#' @param x state vector
#' @param age_cat vector of age categories
#' @param Nvec population distribution (as counts)
#' @examples
#' pp <- read_params("PHAC_testify.csv")
#' ss <- make_state(params=pp)
#' ss2 <- expand_state_age(ss)
#' @export
expand_state_age <- function(x, age_cat=mk_agecats(),
                             Nvec = NULL) {

    ## if no population is provided, assume a uniform distribution
    if(is.null(Nvec)){
      Nvec <- mk_Nvec(age_cat, Ntot = sum(x))
    }

    ## check that number of age groups matches
    if(length(age_cat) != length(Nvec)){
      stop("population distribution must have same length as the number of age categories")
    }

    ## check that state vector and population distribution
    ## sum to the same value
    if(sum(x) != sum(Nvec)){
      stop("state and population distributions must have same sum")
    }

    ## make population distribution
    Ndist <- Nvec/sum(Nvec)

    ## expand state labels with ages
    new_names <- expand_names(names(x), age_cat)

    ## distribute total state counts over subcategories in a way that respects
    ## the overall population distribution
    x <- (
      ## get all states but S
      x[!grepl("^S", names(x))]
      ## distribute counts across all states but S
      ## guided by population distribution
        %>% map_dfr(distribute_counts, dist = Ndist)
      ## get total counts by age over all states but S
        %>% mutate(total = rowSums(across(everything())))
      ## calculate S by age as the age-specific population - count over all
      ## other states
        %>% mutate(S = Nvec - total)
      ## drop total col and move S to far left of the table (so state columns
      ## are ordered as we usually order them)
        %>% select(-total)
        %>% dplyr::relocate(S)
      )

    ## flatten to vector, AND ensure counts
    ## are first sorted by age category, then by state category
    x <- as.vector(t(x))

    ## add sensible names and age cat attribute to output
    names(x) <- new_names
    attr(x, "age_cat") <- age_cat

    ## add class
    class(x) <- "state_pansim"

    return(x)
}

##' collapse age/testing/etc. categories (only; don't do other condensation)
##' @param x a state vector
##' @param levels levels/sort order
##' @export
## wrap into condense() ?
condense_age <- function(x,levels=unique(epi_cat)) {
  ## FIXME: should work on data frames too ...
  epi_cat <- gsub("_.*$","",names(x))
  epi_cat <- factor(epi_cat,levels=levels)
  ret <- vapply(split(x,epi_cat),sum,numeric(1))

  ## re-add class (lost after vapply)
  class(ret) <- "state_pansim"

  return(ret)
}

#' collapse (non-accumulator) states into subcategories (ages)
#' @param x a state vector or data frame (each row is a different time point)
#' @export
condense_state <- function(x){

  ## define non-accumulator state names
  states <- c("S", "E",
              "I", "H", "ICU",
              "R", "D")
  states_regex <- paste0("^", states)

  ## if x is a vector (e.g. a state vec), convert it to a (wide) df
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
    ## separate state from subcategory
    %>% separate(var,
                 into = c("state", "subcat"),
                 sep = "_",
                 extra = "merge")
    ## turn subcat column into factor to keep original ordering
    %>% mutate(subcat = as.factor(subcat))
    ## aggregate value by timestep and subcategory
    %>% group_by(obs_number, subcat)
    %>% summarise(value = sum(value), .groups = "drop")
    ## put data back into wide format
    %>% pivot_wider(names_from = "subcat",
                    values_from = "value")
    ## drop observation number
    %>% select(-obs_number)
  )  -> x

  return(x)
}

## POPULATION DISTRIBUTIONS

#' generate a population distribution (as a vector of counts)
#' @param age_cat age categories (made with `mk_agecats()`)
#' @param Ntot total population size
#' @param dist (character) either "unif" for uniform or "rand" for a random population distribution
#' @param names (logical) return a named vector? (names are given age categories)
#' @export
mk_Nvec <- function(age_cat = mk_agecats(),
                    Ntot = 1e6,
                    dist = "unif",
                    names = FALSE){

  if(!(dist %in% c("unif", "rand"))) stop("dist must be either 'unif' or 'rand'")

  n <- length(age_cat)

  ## uniform distribution
  if(dist == "unif"){
    Ndist <- rep(1/n, n)
  }

  ## random distribution
  if(dist == "rand"){
    Ndraw <- runif(n)
    Ndist <- Ndraw/sum(Ndraw) ## normalize
  }

  ## distribute scale up population distribution to match total population count
  ## in ss
  Nvec <- distribute_counts(total = Ntot, dist = Ndist)

  ## add names?
  if (names){
    names(Nvec) <- age_cat
  }

  return(Nvec)
}

## CONTACT MATRICES

#' generate a contact matrix (where each row is a probability distribution)
#' @param age_cat age categories (made with `mk_agecats()`)
#' @param dist (character) either "unif" for uniform or "rand" for a random distributions
#' @export
mk_Cmat <- function(age_cat = mk_agecats(),
                    dist = "unif"){

  if(!(dist %in% c("unif", "rand"))) stop("dist must be either 'unif' or 'rand'")

  ## set up matrix
  n <- length(age_cat)

  if(dist == "unif"){
    Cmat <- matrix(1/n, nrow=n, ncol=n, dimnames=list(age_cat, age_cat))
  }

  if(dist == "rand"){
    ## preallocate memory
    Cmat <- matrix(nrow = n, ncol = n, dimnames = list(age_cat, age_cat))

    ## fill matrix with a random distribution in each row)
    for (i in 1:n){
      row <- runif(n)
      row <- row/sum(row)

      Cmat[i,] <- row
    }
  }

  return(Cmat)
}

#' Make contact matrix using Mistry et al. approach
#'
#' @param weights named list containing setting-specific weights (in units of average contacts in the given setting per individual of age i with individuals of age j per day)
#' @param province province for which to construct the matrix (if NULL, make a Canada-wide contact matrix)
#'
#' @return matrix of average contacts between individuals of ages i and j per individual of age i per day
#' @export
#'
#' @examples mk_mistry_cmat()
mk_mistry_cmat <- function(weights =
                             c(household = 4.11,
                               school = 11.41,
                               work = 8.07,
                               community = 2.79),
                           province = "Ontario"){

  ## check that weights were specified correctly
  if(!all.equal(sort(names(weights)),
                c("community", "household", "school", "work"))){
    stop("weights vector must be named with names 'household', 'school', 'work', 'community'")
  }

  ## preallocate memory for the output
  cmat <- matrix(rep(0, 85*85), nrow = 85)

  ## combine setting-specific frequency matrices through a
  ## linear combination with the specified weights (in units
  ## of avg number of setting-specific contacts per
  ## individual of age i per unit time) to generate an
  ## overall contact matrix (in units of avg number of
  ## contacts per individual of age i per day)

  ## set up filename prefix/suffix for each
  ## setting-specific matrix
  if(is.null(province)){
    filename_prefix <- paste0("Canada_country_level_F_")
  } else {
    filename_prefix <- paste0("Canada_subnational_", province,
                              "_F_")
  }
  filename_suffix <- "_setting_85.csv"

  settings <- c("household", "school", "work", "community")
  for (set in settings){
    filename <- system.file("params", "mistry-cmats",
                            paste0(filename_prefix,
                                   set,
                                   filename_suffix),
                            package = "McMasterPandemic")
    ## load setting-specific matrix
    set_mat <- readr::read_csv(filename,
                               col_names = FALSE,
                               col_types = cols(
                                 .default = col_double()
                               )) %>%
      as.matrix()

    ## update overall contact matrix by adding a weighted
    ## version of the current setting-specific frequency
    ## matrix
    cmat <- cmat + weights[set]*set_mat
  }

  ## update row and colnames of cmat with age categories
  age_cats <- mk_agecats(min = 0, max = 84, da = 1)
  rownames(cmat) <- age_cats
  colnames(cmat) <- age_cats

  return(cmat)
}

#' expand parameter list to include age structure
#'
#' @param pp parameter list (e.g. read in with `read_params()`)
#' @param age_cat vector of age categories
#' @param Cmat contact matrix; default is uniform
#' @param Nvec population distribution (as counts); default is uniform
#' @examples
#' pp <- read_params("PHAC_testify.csv")
#' ppa <- expand_params_age(pp)
#' @export
expand_params_age <- function(pp,
                              age_cat = mk_agecats(),
                              Cmat = mk_Cmat(),
                              Nvec = mk_Nvec()){
  ## convert to list
  pp <- as.list(pp)
  ## update pop
  pp[["N"]] <- Nvec
  ## update Cmat
  pp <- c(pp, list(Cmat = Cmat))

  ## add age cats as attribute
  attr(pp, "age_cat") <- age_cat

  ## add class
  class(pp) <- "params_pansim"

  return(pp)
}

#' @rdname expand_state_age
#'
#' @param ratemat rate matrix
#' @param params parameter vector
#' @examples
#' pp <- read_params("PHAC_testify.csv")
#' state <- make_state(params=pp)
#' M <- make_ratemat(state,pp, sparse=TRUE)
#' Ma <- ageify(M, pp)
#' library(Matrix)
#' Matrix::image(Ma)
#' Mta <- ageify(testify(M,pp),pp)
#' Matrix::image(Mta)
#' @export
ageify <- function(ratemat, params, age_cat=mk_agecats()) {
    m <- Matrix::kronecker(diag(length(age_cat)), ratemat)
    new_names <- unlist(lapply(rownames(m), paste, age_cat, sep="_"))
    dimnames(m) <- list(new_names, new_names)
    return(m)
}

## FIXME: age-dependent params??

if (FALSE) {
    devtools::load_all()
    pp <- read_params("PHAC_testify.csv")
    ss <- make_state(params=pp)
    ss2 <- expand_state_age(ss)
    ## hack so we have an infective
    ss2["Im_11-20"] <- 1
    ss2["E_91+"] <- 0
    tot_I <- function(x) sum(x[grep("^I[a-z]",names(x))])
    tot_I(ss)
    sum(ss2)
    tot_I(ss2)
    condense.pansim(data.frame(date=NA,rbind(ss2)),add_reports=FALSE)
    M <- make_ratemat(ss2, pp, sparse=TRUE)
    show_ratemat(M)
    aa <- mk_agecats()
    ## compound symmetric example
    Cmat <- matrix(0.1, nrow=length(aa), ncol=length(aa), dimnames=list(aa,aa))
    diag(Cmat) <- 1
    ifun <- function(M) {
        Matrix::image(Matrix(M),scales=list(y=list(at=seq(nrow(M)),labels=rownames(M)),
                                            x=list(at=seq(ncol(M)),labels=colnames(M), rot=90)))
    }
    ppa <- c(as.list(pp),list(Cmat=Cmat))
    b1 <- make_betavec(ss2, ppa, full=FALSE)
    ifun(b1)
    b2 <- Matrix(make_betavec(ss2, ppa, full=TRUE))
    ifun(b2)
    M <- make_ratemat(ss2, ppa, sparse=TRUE)
    show_ratemat(M)
    M %*% ss2
}

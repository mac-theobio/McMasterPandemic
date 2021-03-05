#' construct age categories
#' @param min minimum age
#' @param max maximum age
#' @param da age bin width
#' @export
mk_agecats <- function(min=1,max=100,da=10) {
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
    return(ret)
}

## x_y, with x varying faster
expand_names <- function(x,y,sep="_") {
    unlist(lapply(y, function(a) paste(x, a, sep=sep)))
}

#' distribute counts given a desired distribution
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

#' expand state vector and rate matrix by age classes and population distribution
#'
#' epidemiological state varies fast, age category varies slowly
#' @param x state vector
#' @param age_cat vector of age categories
#' @param N_dist distribution of population over ages (sums to 1)
#' @examples
#' pp <- read_params("PHAC_testify.csv")
#' ss <- make_state(params=pp)
#' ss2 <- expand_stateval_age(ss)
#' @export
expand_stateval_age <- function(x, age_cat=mk_agecats(),
                                N_dist = NULL) {
  
    ## if no population is provided, assume a uniform distribution
    if(is.null(N_dist)){
      tot_N <- sum(x)
      n_agecats <- length(age_cat)
      N_dist <- rep(tot_N/n_agecats, n_agecats)/tot_N
    }
    
    ## check that population distribution vector sums to 1
    if(sum(N_dist) != 1){
      stop("the population distribution (N_dist) must sum to 1")
    }
  
    ## check that number of age groups matches 
    if(length(age_cat) != length(N_dist)){
      stop("state and population vectors must have same length")
    }
    
    ## expand state labels with ages
    new_names <- expand_names(names(x), age_cat)
    
    ## split total count for each class over age categories
    ## based on population distribution
    ## (map_dfr -> t -> as vector business is to ensure counts
    ## are returned sorted by age category, then by state category)
    x <- as.vector(t(map_dfr(x, distribute_counts, dist = N_dist)))
    names(x) <- new_names

    attr(x, "age_cat") <- age_cat
    return(x)
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

## FIXME: carry age categories as attribute of stateval?
## assign class state_pansim?

#' @rdname expand_stateval_age
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
    ss2 <- expand_stateval_age(ss)
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

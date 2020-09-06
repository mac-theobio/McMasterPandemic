mk_agecats <- function(min=1,max=100,da=10) {
    s1 <- seq(min,max,by=da)
    c(sprintf("%d-%d",s1[-length(s1)],(s1[-1]-1)),
      paste0(s1[length(s1)],"+"))
}

#' expand state vector and rate matrix by age classes
#' 
#' @examples
#' pp <- read_params("PHAC_testify.csv")
#' ss <- make_state(params=pp)
#' ss2 <- expand_stateval_age(ss)
expand_stateval_age <- function(x, age_cat=mk_agecats()) {
    new_names <- unlist(lapply(names(x), paste, age_cat, sep="_")) ## or outer() ?
    n_expand <- length(age_cat)
    new_states <- smart_round(rep(x, each=n_expand)/n_expand)
    names(new_states) <- new_names
    return(new_states)
}

## FIXME: carry age categories as attribute of stateval?
## assign class state_pansim?

#' @rdname expand_stateval_age
#' 
#' @examples
#' pp <- read_params("PHAC_testify.csv")
#' state <- make_state(params=pp)
#' M <- make_ratemat(state,pp, sparse=TRUE)
#' Ma <- ageify(M, pp)
#' library(Matrix)
#' Matrix::image(Ma)
#' Mta <- ageify(testify(M,pp),pp))
#' Matrix::image(Mta)
ageify <- function(ratemat, params, age_cat=mk_agecats()) {
    m <- Matrix::kronecker(diag(length(age_cat)), ratemat)
    new_names <- unlist(lapply(rownames(m), paste, age_cat, sep="_"))
    dimnames(m) <- list(new_names, new_names)
    return(m)
}
## FIXME: age-dependent params??
    

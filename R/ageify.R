#' construct age categories
#' @param min minimum age
#' @param max maximum age
#' @param da age bin width
#' @export
mk_agecats <- function(min=1,max=100,da=10) {
    s1 <- seq(min,max,by=da)
    c(sprintf("%d-%d",s1[-length(s1)],(s1[-1]-1)),
      paste0(s1[length(s1)],"+"))
}

## x_y, with x varying faster
expand_names <- function(x,y,sep="_") {
    unlist(lapply(y, function(a) paste(x, a, sep=sep)))
}

#' expand state vector and rate matrix by age classes
#'
#' epidemiological state varies fast, age category varies slowly
#' @param x state vector
#' @param age_cat vector of age categories
#' @examples
#' pp <- read_params("PHAC_testify.csv")
#' ss <- make_state(params=pp)
#' ss2 <- expand_stateval_age(ss)
#' @export
expand_stateval_age <- function(x, age_cat=mk_agecats()) {
    new_names <- expand_names(names(x), age_cat)
    n_expand <- length(age_cat)
    new_states <- rep(x, n_expand)/n_expand
    names(new_states) <- new_names
    ## round Susceptible and non-Susceptible compartments, maintaining sum *separately*
    ## so we don't lose all the non-susc ...
    S_pos <- grep("^S", names(new_states))
    new_states[S_pos] <- smart_round(new_states[S_pos])
    new_states[-S_pos] <- smart_round(new_states[-S_pos])
    return(new_states)
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

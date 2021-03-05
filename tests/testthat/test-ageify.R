library(McMasterPandemic)
## devtools::load_all()
library(testthat)

## apparently context() is depreciated/superseded
## tests are now be auto-named with the filename
## ?testthat::context
context("ageify")

## set up testing environment ##

## base params and states 
pp <- update(read_params("PHAC_testify.csv"), testing_intensity=0)
ss <- make_state(params=pp)

## set age categories
age_cat <- mk_agecats()
## draw random population distribution
N_draw <- rpois(length(age_cat), lambda = 10)
N_dist <- N_draw/sum(N_draw) ## normalize as a distribution
## distribute scale up population distribution to match total population count
## in ss
N_vec <- distribute_counts(total = sum(ss), dist = N_dist)
## generate state vec
ss2 <- expand_stateval_age(ss, age_cat = age_cat, N_vec = N_vec)

## tests ##

test_that("population sizes don't change when adding age structure",
{
    ## total population sizes
    expect_equal(sum(ss2), sum(ss))
    ## population size of each state
    state_regex <- paste0("^", names(ss))
    counts_by_state <- unlist(map(state_regex, total_by_cat, x = ss2))
    expect_equal(counts_by_state,
                 as.vector(ss))
    ## population size of each age class
    age_cat <- attr(ss2, "age_cat")
    age_cat <- sub("\\+", "\\\\+", age_cat) ## escape + in age group for regex
    age_cat_regex <- paste0(age_cat,"$")
    counts_by_age <- unlist(map(age_cat_regex, total_by_cat, x = ss2))
    expect_equal(counts_by_age,
                 N_vec)
})

## not really proper tests yet: FIXME/clean me up!
test_that("generic age stuff", {
    #condense.pansim(data.frame(date=NA,rbind(ss2)),add_reports=FALSE)
    M <- make_ratemat(ss2, pp, sparse=TRUE)
    show_ratemat(M)
    ## convert params to list in prep for adding population
    ## size vector and contact matrix update population
    ppa <- as.list(pp)
    ## param to a vector of age-specific populations
    ## (uniform distribution)
    N_vec <- rep(pp[["N"]]/length(age_cat), length(age_cat))
    ppa <- update(ppa, N = N_vec)
    ## compound symmetric example (for a uniformly
    ## distributed population---otherwise Cmat would not be
    ## symmetric since each row should be divided by the
    ## size of the susceptible population)
    Cmat <- matrix(0.1, nrow=length(age_cat), ncol=length(age_cat),
                   dimnames=list(age_cat,age_cat))
    diag(Cmat) <- 1
    ppa <- c(ppa,list(Cmat=Cmat))
    ifun <- function(M) {
        Matrix::image(Matrix(M),scales=list(y=list(at=seq(nrow(M)),labels=rownames(M)),
                                            x=list(at=seq(ncol(M)),labels=colnames(M), rot=90)))
    }
    b1 <- make_betavec(ss2, ppa, full=FALSE)
    expect_equal(dim(b1), c(10,40))
    ifun(b1)
    b2 <- Matrix(make_betavec(ss2, ppa, full=TRUE))
    ifun(b2)
    expect_equal(dim(b2), c(10,130))
    M <- make_ratemat(ss2, ppa, sparse=TRUE)
    expect_equal(dim(M),c(130,130))
    show_ratemat(M)
    M %*% ss2
    rr <- run_sim_range(ppa, ss2, nt=1000)
    ## suppress warnings about values <0 (??)
    suppressWarnings(matplot(rr[,1],rr[,-1],lty=1,type="l",log="y"))
    rr2 <- run_sim(ppa, ss2,end_date="2022-Nov-1",condense=FALSE)
    plot(rr2,log=TRUE)+theme(legend.position="none")
})

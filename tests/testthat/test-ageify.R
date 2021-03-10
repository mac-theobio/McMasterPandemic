library(McMasterPandemic)
## devtools::load_all()
library(testthat)

## apparently context() is depreciated/superseded
## tests are now be auto-named with the filename
## ?testthat::context
context("ageify")

## SETUP ##

## base params and state
params <- update(read_params("PHAC_testify.csv"), testing_intensity=0)
state <- make_state(params=params)

## generate state vecs ##

## uniform population distribution
Nvec_u <- mk_Nvec(Ntot = sum(state))
state_u <- expand_state_age(state)

## random population distribution
Nvec_r <- mk_Nvec(Ntot = sum(state),
                     dist = "rand")
state_r <- expand_state_age(state, Nvec = Nvec_r)

## generate param sets ##

Cmat_r <- mk_Cmat(dist = "rand")

## params with unif pop
## + unif cmat
params_uu <- expand_params_age(params)
## + random cmat
params_ur <- expand_params_age(params,
                               Cmat = Cmat_r)

## params with random pop
## + unif cmat
params_ru <- expand_params_age(params,
                               Nvec = Nvec_r)
## + random cmat
params_rr <- expand_params_age(params,
                               Nvec = Nvec_r,
                               Cmat = Cmat_r)

## generate sims ###
end_date <- "2021-02-15"

## sims with unif pop
## + unif cmat
res_uu <- run_sim(params_uu, state_u, end_date = end_date, condense = FALSE)
## + rand cmat
res_ur <- run_sim(params_ur, state_u, end_date = end_date, condense = FALSE)

## sims with random pop
## + unif cmat
res_ru <- run_sim(params_ru, state_r, end_date = end_date, condense = FALSE)
## + rand cmat
res_rr <- run_sim(params_rr, state_r, end_date = end_date, condense = FALSE)

## TESTS

## population count tests
    
test_that("population distributions are properly initialized", {
    ## uniform population
    expect_equal(sum(mk_Nvec(Ntot = 1e6)), 1e6)
    ## random population
    expect_equal(sum(mk_Nvec(Ntot = 1e6, dist = "rand")), 1e6)
})

test_that("initial state population sizes don't change after adding age structure",
{
    ## total population sizes
    expect_equal(sum(state_u), sum(state))
    expect_equal(sum(state_r), sum(state))

    ## population size of each state
    expect_equal(condense_age(state_u),
                 state)
    expect_equal(condense_age(state_r),
                 state)
    
    ## population size of each age class
    ## uniform population
    expect_equal(as.numeric(condense_state(state_u)),
                 Nvec_u)
    ## random population
    expect_equal(as.numeric(condense_state(state_r)),
                 Nvec_r)
})

## helper function to check that population remains constant across age groups
## and time steps in a simulation
check_const_pop <- function(res, params){
    ## condense states
    res_pops <- (condense_state(res)
                 ## convert rows to a single list-col containing the age-specific
                 ## population distribution at each time step
                 %>% transmute(dist = transpose(select(.,everything()))))
    
    ## check every row of the sim result data frame against the pop distribution
    ## (doing ifelse here because i don't
    ## know how to get all.equal to return FALSE instead of mean relative diff)
    check_rows <- map_lgl(res_pops$dist,
                          ~ isTRUE(all.equal(unname(unlist(.)), params[["N"]])))
    ## check that all rows passed the test 
    check_all <- all(check_rows)
    return(check_all)
}

test_that("age-specific population doesn't change over the course of a simulation",
{
    ## unif pop, unif Cmat
    expect_true(check_const_pop(res_uu, params_uu))
    
    ## unif pop, rand Cmat
    expect_true(check_const_pop(res_ur, params_ur))
    
    ## rand pop, unif Cmat
    expect_true(check_const_pop(res_ru, params_ru))
    
    ## rand pop, rand Cmat
    expect_true(check_const_pop(res_rr, params_rr))
})

## beta tests

test_that("age-structured beta0 has correct dimensions", {
    expect_identical(dim(make_betavec(state_u, params_uu)),
                     as.integer(c(1, 13))*length(mk_agecats()))
})

## simulation tests

test_that("homogeneous case of age-structured model reduces to base (non-ageified) model (comparing simulations)", {
    ## condense homogeneous simulation
    res_uu_cond <- condense.pansim(res_uu)
    ## base sim (no age-structure, same params)
    res_hom <- run_sim(params, state, end_date = end_date)
    ## check subset of state variables (some special cols in sim not set up for
    ## ageify case, like foi)
    expect_equal((res_uu_cond %>% select(S:D)),
                 (res_hom %>% select(S:D)))
})

## not really proper tests yet: FIXME/clean me up!
# test_that("generic age stuff", {
#     b1 <- make_betavec(ss2, ppa, full=FALSE)
#     expect_equal(dim(b1), c(10,40))
#     ifun(b1)
#     b2 <- Matrix(make_betavec(ss2, ppa, full=TRUE))
#     ifun(b2)
#     expect_equal(dim(b2), c(10,130))
#     M <- make_ratemat(ss2, ppa, sparse=TRUE)
#     expect_equal(dim(M),c(130,130))
#     show_ratemat(M)
#     M %*% ss2
#     rr <- run_sim_range(ppa, ss2, nt=1000)
#     ## suppress warnings about values <0 (??)
#     suppressWarnings(matplot(rr[,1],rr[,-1],lty=1,type="l",log="y"))
#     rr2 <- run_sim(ppa, ss2,end_date="2022-Nov-1",condense=FALSE)
#     plot(rr2,log=TRUE)+theme(legend.position="none")
# })

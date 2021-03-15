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
## set up so there will be a single infectious seed in each age group (for the
## age-structured sims)
state["E"] <- length(mk_agecats())
state["S"] <- 1e6 - state["E"]
state[!grepl("^(S|E)", names(state))] <- rep(0, length(state)-2)

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
Cmat_d <- mk_Cmat(dist = "diag")

## params with unif pop
## + unif cmat
params_uu <- expand_params_age(params)
## + random cmat
params_ur <- expand_params_age(params,
                               Cmat = Cmat_r)
## + diag cmat
params_ud <- expand_params_age(params,
                               Cmat = Cmat_d)

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
## + diag cmat
res_ud <- run_sim(params_ud, state_d, end_date = end_date, condense = FALSE)

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
                     as.integer(c(1, length(state)))*length(mk_agecats()))
})

## ratemat tests

test_that("age-structure rate matrix has correct dimensions", {
    M <- make_ratemat(state_u, params_uu, sparse=TRUE)
    expect_equal(dim(M),
                 c(length(state), length(state))*length(mk_agecats()))
})

## simulation tests

test_that("homogeneous case of age-structured model condenses to base (non-ageified) model (comparing simulations)", {
    ## condense homogeneous simulation
    res_uu_cond <- condense.pansim(res_uu)
    ## base sim (no age-structure, same params)
    res_hom <- run_sim(params, state, end_date = end_date)
    ## check subset of state variables (some special cols in sim not set up for
    ## ageify case, like foi)
    expect_equal((res_uu_cond %>% select(S:D)),
                 (res_hom %>% select(S:D)))
})

## utility function to check that all elements of a vector are equal
## based on:
## https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-numeric-vector
zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
    if (length(x) == 1) return(TRUE)
    if(anyNA(x)) return(FALSE)
    ## if the input is a list, unpack it
    if(is.list(x)) x <- unlist(x)
    x <- range(x) / mean(x)
    return(isTRUE(all.equal(x[1], x[2], tolerance = tol)))
}

## helper function to check for equality across age groups for a given state and
## simulation date (returns a logical vector of length date*states)
check_equality_across_ages <- function(df){
    (df
        ## pivot to be able to easily group observations using substrings in column
        ## name
        %>% select(-foi)
        %>% pivot_longer(-date)
        %>% separate(name, into = c("state", "age_cat"),
                     sep = "_", extra = "merge")
        ## expand age classes back into separate columns,
        ## and then compress into a single list-col,
        ## so that we can iterate over a list at a time to check for no variance
        ## "zero range"
        %>% pivot_wider(names_from = "age_cat")
        %>% mutate(values = transpose(select(.,matches("\\d+"))))
        ## check for equality across age groups
        %>% mutate(values_all_equal = map_lgl(values, zero_range))
        ## pull check vector to use in an expectation
        %>% pull(values_all_equal)
     ) -> values_all_equal

    return(values_all_equal)
}

test_that("homogeneous case of age-structured model yields identical epidemics in age classes that did not seed the epidemic",{
    expect_true(all(check_equality_across_ages(res_uu)))
})

test_that("diagonal contacts yield identical epidemics in each age group (with an infectious seed in each age group)",{
    expect_true(all(check_equality_across_ages(res_ud)))
})

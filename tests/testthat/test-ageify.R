library(McMasterPandemic)
# devtools::load_all()
library(testthat)
library(purrr) ## for transpose()

## apparently context() is depreciated/superseded
## tests are now be auto-named with the filename
## ?testthat::context
context("ageify")

## SETUP ##

## base params and state
age_cat <- mk_agecats()
beta0 <- 2
params <- update(read_params("PHAC_testify.csv")
                 , testing_intensity=0
                 , beta0 = beta0
                 )
state <- make_state(N = params[["N"]],
                    E0 = length(age_cat))

## generate population distributions
Nvec_u <- mk_Nvec(Ntot = params[["N"]], age_cat = age_cat,
                  dist = "unif")
## random population distribution
# Nvec_r <- mk_Nvec(Ntot = sum(state), age_cat = age_cat,
#                   dist = "rand")

## generate state vecs ##
state_u <- expand_state_age(state, age_cat = age_cat,
                            Nvec = Nvec_u)
# state_r <- expand_state_age(state, age_cat = age_cat,
#                             Nvec = Nvec_r)

## generate param sets ##

## params with unif pop
## + unif pmat
params_uu <- expand_params_age(params, age_cat = age_cat,
                               pmat = mk_pmat(age_cat = age_cat,
                                              dist = "unif"),
                               Nvec = Nvec_u)
## + random pmat
# params_ur <- expand_params_age(params, age_cat = age_cat,
#                                pmat = mk_pmat(age_cat = age_cat,
#                                               dist = "rand"),
#                                Nvec = Nvec_u)
## + diag pmat
params_ud <- expand_params_age(params, age_cat = age_cat,
                               pmat = mk_pmat(age_cat = age_cat,
                                              dist = "diag"),
                               Nvec = Nvec_u)

## params with random pop
## + unif pmat
# params_ru <- expand_params_age(params, age_cat = age_cat,
#                                pmat = mk_pmat(age_cat = age_cat,
#                                               dist = "unif",
#                                               Nvec = Nvec_r),
#                                Nvec = Nvec_r)
## + random pmat
# params_rr <- expand_params_age(params, age_cat = age_cat,
#                                pmat = mk_pmat(age_cat = age_cat,
#                                               dist = "rand",
#                                               Nvec = Nvec_r),
#                                Nvec = Nvec_r)

## generate sims ###
end_date <- "2021-02-15"

## sims with unif pop
## + unif pmat
res_uu <- run_sim(params_uu, state_u, end_date = end_date,
                  condense_args = c(keep_all = TRUE)) ## don't condense age groups
## + rand pmat
# res_ur <- run_sim(params_ur, state_u, end_date = end_date, condense = FALSE)
## + diag pmat
res_ud <- run_sim(params_ud, state_u, end_date = end_date,
                  condense_args = c(keep_all = TRUE))

## sims with random pop
## + unif pmat
# res_ru <- run_sim(params_ru, state_r, end_date = end_date, condense = FALSE)
## + rand pmat
# res_rr <- run_sim(params_rr, state_r, end_date = end_date, condense = FALSE)

## TESTS

## preprocessing steps
#########################

## population count

test_that("population distributions are properly initialized", {
    ## uniform population
    expect_equal(sum(mk_Nvec(Ntot = 1e6)), 1e6)
    ## random population
    # expect_equal(sum(mk_Nvec(Ntot = 1e6, dist = "rand")), 1e6)
})

test_that("initial state population sizes don't change after adding age structure",
{
    ## total population sizes
    expect_equal(sum(state_u), sum(state))
    # expect_equal(sum(state_r), sum(state))

    ## population size of each state
    # expect_equal(condense_age(state_u),
    #              state)
    # expect_equal(condense_age(state_r),
    #              state)

    ## population size of each age class
    ## uniform population
    expect_equal(as.numeric(condense_state(state_u)),
                 Nvec_u)
    ## random population
    # expect_equal(as.numeric(condense_state(state_r)),
    #              Nvec_r)
})

## expand_params_age
test_that("ageified parameters are correctly initialized",{
    ## defaults
    expect_identical(expand_params_age(params, age_cat = age_cat)$beta0, beta0)
    ## update beta0 manually
    test_vec <- 1:length(age_cat)
    names(test_vec) <- age_cat
    expect_identical(expand_params_age(
        params, age_cat = age_cat, beta0 = test_vec,
        balance_warning = FALSE)$beta0,
        test_vec)
    ## specify transmissibility and contact_rate_age, calculate beta0
    expect_identical(expand_params_age(
        params, age_cat = age_cat,
        transmissibility = 0.5,
        contact_rate_age = test_vec,
        balance_warning = FALSE)$beta0,
        0.5*test_vec)
})

## pmat
test_that("Mistry et al. contact matrix get aggregated with the correct dims",{
    age_cat <- mk_agecats(min = 0, max = 84, da  = 10)
    mistry_pmat <- expand_params_mistry(params, age_cat = age_cat)$pmat
    expect_equal(dim(mistry_pmat),
                 c(length(age_cat), length(age_cat)))
})

## beta
test_that("age-structured beta0 has correct dimensions", {
    expect_identical(dim(make_beta(state_u, params_uu)),
                     as.integer(c(1, length(state)))*length(mk_agecats()))
})

## ratemat
test_that("age-structure rate matrix has correct dimensions", {
    M <- make_ratemat(state_u, params_uu, sparse=TRUE)
    expect_equal(dim(M),
                 c(length(state), length(state))*length(mk_agecats()))
})

## simulation tests
######################

## helper function to check that population remains constant across age groups
## and time steps in a simulation
check_const_pop <- function(res, params){
    ## condense states
    res_pops <- (condense_state(res)
                 ## convert rows to a single list-col containing the
                 ## age-specific population distribution at each time step
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
              ## unif pop, unif pmat
              expect_true(check_const_pop(res_uu, params_uu))

              ## unif pop, rand pmat
              # expect_true(check_const_pop(res_ur, params_ur))

              ## rand pop, unif pmat
              # expect_true(check_const_pop(res_ru, params_ru))

              ## rand pop, rand pmat
              # expect_true(check_const_pop(res_rr, params_rr))
          })

test_that("homogeneous case of age-structured model condenses to base (non-ageified) model (comparing simulations)", {
    ## condense homogeneous simulation, get rid of any cols that still have age
    ## groups post condense
    res_uu_cond <- condense(res_uu)
    ## base sim (no age-structure, same params)
    res_hom <- run_sim(params, state, end_date = end_date)
    ## check subset of state variables (some special cols in sim not set up for
    ## ageify case, like foi)
    expect_equal((res_uu_cond %>% select(-starts_with("foi"))),
                 (res_hom %>% select(-foi)))
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
        %>% select(!starts_with("foi"))
        %>% pivot_longer(-date)
        %>% separate(name, into = c("state", "age_cat"),
                     sep = "_", extra = "merge")
        ## expand age classes back into separate columns,
        ## and then compress into a single list-col,
        ## so that we can iterate over a list at a time to check for no variance
        ## "zero range"
        %>% pivot_wider(names_from = "age_cat")
        %>% mutate(values = transpose(select(.,matches("\\d+"))))
        ## remove observations where all age-specific values are NA
        %>% filter(map(map(values, is.na), sum)==0)
        ## check for equality across age groups
        %>% mutate(values_all_equal = map_lgl(values, zero_range))
        ## pull check vector to use in an expectation
        %>% pull(values_all_equal)
     ) -> values_all_equal

    return(values_all_equal)
}

drop_agg_reports <- function(df){
    df <- (df
        %>% select(!c(incidence, report, cumRep)))

    return(df)
}

test_that("homogeneous case of age-structured model yields identical epidemics in age classes that did not seed the epidemic",{
    ## remove aggregate reporting columns
    expect_true(all(check_equality_across_ages(drop_agg_reports(res_uu))))
})

test_that("diagonal contacts yield identical epidemics in each age group (with an infectious seed in each age group)",{
    expect_true(all(check_equality_across_ages(drop_agg_reports(res_ud))))
})

test_that("Mistry contact parameters are properly initialized", {
    ## using original age groups
    Nvec <- mk_mistry_Nvec(province = "Alberta")
    params_mistry <- expand_params_mistry(params,
                                          province = "Alberta")
    ## check population distribution
    expect_equal(params_mistry$N, Nvec)
    ## check that we recover a symmetric contact abundance matrix from
    ## mistry frequency matrices (as stored)
    expect_true(isSymmetric(params_mistry$mistry_fmats$community*params_mistry$N))
    ## check pmat rows sum to 1
    expect_equal(unname(rowSums(params_mistry$pmat)),
                 rep(1, length(attr(params_mistry, "age_cat"))))

    ## using aggregated age groups
    age_cat <- mk_agecats(min = 0, max = 80, da = 10)
    Nvec <- mk_mistry_Nvec(province = "Quebec", age_cat = age_cat)
    params_mistry <- expand_params_mistry(params,
                                          province = "Quebec",
                                          age_cat = age_cat)

    ## check population distribution
    expect_equal(params_mistry$N, Nvec)
    ## check that we recover a symmetric contact abundance matrix from
    ## mistry frequency matrices (as stored)
    expect_true(isSymmetric(params_mistry$mistry_fmats$community*params_mistry$N))

})

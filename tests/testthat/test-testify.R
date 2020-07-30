library(McMasterPandemic)
library(testthat)

context("testify")

pp <- read_params("PHAC_testify.csv")

## Making states and expanding states
state <- make_state(params=pp)
state_testified <- expand_stateval(state)

test_that("testified states make sense", {
    expect_equal(length(state_testified), 4*length(state)+2)
    expect_equal(sort(unique(gsub("_.*$","",names(state_testified)))),
                 c(sort(c("N","P",names(state)))))
})

## Making beta_vec wtr states
beta_vec <- make_betavec(state,pp)
beta_vec_testified <- make_betavec(state_testified,pp)

test_that("testified betas make sense", {
    expect_equal(names(beta_vec_testified), names(state_testified))
    expect_equal(sort(unique(beta_vec_testified)),
                 sort(unique(beta_vec)))
})

test_that("catch state/beta mismatch", {
    expect_error(update_foi(state,pp,beta_vec_testified), "not the same")
})

## Making ratemat
ratemat <- make_ratemat(state,pp)
ratemat_testified <- testify(ratemat,pp) 
test_that("ratemat makes sense", {
    expect_equal(ncol(ratemat_testified), length(state_testified))
})

## Updating FOI
test_that("FOI doesn't change", {
    expect_equal(update_foi(state,pp,beta_vec),
                 update_foi(state_testified,pp,beta_vec_testified))
})

test_that("condensation is OK", {
    sim0_testified_condensed <- run_sim(params = pp,
                                        ratemat_args = list(testify=TRUE))
    expect_equal(names(sim0_testified_condensed),
                 c("date", "S", "E", "I", "H", "hosp",
                   "ICU", "R", "death", "foi", 
                   "incidence", "report", "cumRep", "D"))
})

library(testthat)
library(McMasterPandemic)

params <- read_params("ICU1.csv")
state <- make_state(params[["N"]],E0=params[["E0"]], use_eigvec=FALSE)
M <- make_ratemat(state, params)

test_that("basic show_ratemat works", {
  ## FIXME: would like to do better testing but how do we test graphics ... ?
  ## expect_snapshot?
  expect_is(show_ratemat(M), "trellis")
  expect_is(show_ratemat(M, block_size=c(3,5), block_col=c(2,4)), "trellis")
  expect_is(show_ratemat(M, block_size=NA), "trellis")
})

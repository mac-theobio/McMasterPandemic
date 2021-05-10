library(McMasterPandemic)
library(testthat)

local_edition(3)
test_that('Identical .rda', {
load("C:/Users/somat/Documents/GitHub/McMasterPandemic/refactor/ontario-rda/ONcalib_2020Jun01.rda")
expect_snapshot(ont_cal1)

load("C:/Users/somat/Documents/GitHub/McMasterPandemic/data/ONcalib.rda")
expect_snapshot(ont_cal1)

})

test_that('Identical 2brks .rda', {
  load("C:/Users/somat/Documents/GitHub/McMasterPandemic/refactor/ontario-rda/ONcalib_2brks_2020Jun01.rda")
  expect_snapshot(ont_cal1)

  load("C:/Users/somat/Documents/GitHub/McMasterPandemic/data/ONcalib_2brks.rda")
  expect_snapshot(ont_cal1)

})

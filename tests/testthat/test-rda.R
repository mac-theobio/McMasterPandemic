library(McMasterPandemic)
library(testthat)

local_edition(3)
test_that('Identical .rda', {
load(system.file("testdata","ONcalib_2020Jun01.rda",package="McMasterPandemic"))
expect_snapshot(ont_cal1)

load(system.file("data","ONcalib.rda",package="McMasterPandemic"))
expect_snapshot(ont_cal1)

})

test_that('Identical 2brks .rda', {
  load(system.file("testdata","ONcalib_2brks_2020Jun01.rda",package="McMasterPandemic"))
  expect_snapshot(ont_cal1)


  load(system.file("data","ONcalib_2brks.rda",package="McMasterPandemic"))
  expect_snapshot(ont_cal1)

})

library(testthat)
library(McMasterPandemic)
library(ggplot2)

params <- read_params("ICU1.csv")
s <- run_sim(params, start_date="1-Mar-2020", end_date="1-May-2020")

context("aggregation")

test_that("Jacobian/r/etc", {
    J <- make_jac(params)
    expect_equal(unname(colSums(J)), rep(0,nrow(J)))
    expect_equal(get_r(params,"kernel"), get_r(params, "expsim"), tolerance=1e-5)
    if (FALSE) {
        expect_equal(get_r(params,"expsim"), get_r(params, "analytical"), tolerance=1e-5)
        ## still don't know why dominant eigenvalue is wrong! what else to look at to diagnose it?
        ## figure out proper comparisons to rate matrix (which is what actually drives the sim)
        ## M <- make_ratemat(make_state(N=params[["N"]],E=1e-3),params)
    }
})

test_that("basic aggregation", {
    c1 <- condense(s)
    expect_error(condense(s,junk=TRUE), "unknown arguments")
    expect_equal(dim(c1),c(62,11))
    expect_equal(names(c1),
                 c("date", "S", "E", "I", "H", "ICU", "R", "death", "foi",
                   "incidence", "report"))
    first <<- dplyr::first
    expect_error(aggregate(s,junk=TRUE), "unknown arguments")
    a1 <- aggregate(condense(s), start="12-Feb-2020", period="7 days",
                    FUN=list(mean=c("H","ICU","I"),
                             sum=c("report","death")))
    expect_equal(dim(a1), c(10,6))
})

test_that("trans_labels", {
    vv <- c("H", "ICU", "Ventilator", "report", "newTests", "death")
    expect_equal(unique((tv <- trans_state_vars(ont_all))$var),
                 vv)
    ## idempotent ..
    expect_equal(unique(trans_state_vars(tv)$var), vv)
})

test_that("fit methods", {
    suppressWarnings(plot(ont_cal1))
    suppressWarnings(plot(ont_cal1,data=trans_state_vars(ont_all)))
    suppressWarnings(plot(ont_cal_2brks,data=trans_state_vars(ont_all)))
})

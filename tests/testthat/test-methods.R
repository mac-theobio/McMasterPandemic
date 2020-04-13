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
    expect_equal(dim(c1),c(62,10))
    expect_equal(names(c1),
                 c("date", "S", "E", "I", "H", "ICU", "R", "D",
                   "incidence", "report"))
    first <<- dplyr::first
    a1 <- aggregate(condense(s), start="12-Feb-2020", period="7 days",
                    FUN=list(mean=c("H","ICU","I"),
                             first=c("D"),sum=c("report")))

    expect_equal(dim(a1), c(10,6))
})

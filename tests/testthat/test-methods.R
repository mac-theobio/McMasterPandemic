library(testthat)
library(McMasterPandemic)
library(ggplot2)

params <- read_params("ICU1.csv")
s <- run_sim(params, start_date="1-Mar-2020", end_date="1-May-2020")

context("aggregation")

test_that("Jacobian/r/etc", {
    J <- make_jac(params)
    expect_equal(unname(colSums(J)), rep(0,nrow(J)))
    expect_equal(get_r(params,"kernel"), get_r(params, "expsim"), tolerance=2e-3)  ## FIXME: should be closer?
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
    expect_equal(dim(c1),c(62,15))
    expect_equal(names(c1),
                 c("date", "S", "E", "I", "H", "ICU", "R", "hosp", "X",
                   "death", "D", "foi", "incidence", "report", "cumRep"))
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
    expect_is(suppressWarnings(plot(ont_cal1)), "ggplot")
    expect_is(suppressWarnings(plot(ont_cal1,data=trans_state_vars(ont_all))), "ggplot")
    predict(ont_cal_2brks)
    expect_is(suppressWarnings(plot(ont_cal_2brks,data=trans_state_vars(ont_all))), "ggplot")
})

test_that("predict", {
    ## n.b. this will change when we update! add "D"? (hacked for now)
    expect_equal(setdiff(unique(predict(ont_cal1)$var),"D"),
                 c("H","ICU","hosp","death","incidence","report","cumRep"))
    pp0 <- predict(ont_cal1,keep_vars="all",sim_args=list(condense=FALSE))
    pp0_v <- unique(pp0$var)
    expect_equal(length(pp0_v),14L)
    pp1 <- predict(ont_cal1,keep_vars="all")
    pp1_v <- unique(pp1$var)
    expect_equal(length(pp1_v),14L)
    pp2 <- predict(ont_cal1,stoch=c(proc=TRUE,obs=TRUE),
                   stoch_start=c(proc="2020-04-10",obs="2020-01-30"),
                   new_params=c(proc_disp=5,obs_disp=100))
    expect_is(pp2,"predict_pansim")
    plot(pp2)
    pred_Rt <- predict(ont_cal1, keep_vars="Rt")
    ## don't want to test exact values, because ont_cal1 will change
    ## with new data: should fit to simulated data and test that ...
    suppressWarnings(pp2 <- predict(ont_cal1, ensemble=TRUE,
                                    .progress="none",
                                    imp_wts=TRUE, nsim=10))
    expect_equal(dim(pp2),c(861,6))
    set.seed(101)
    suppressWarnings(pp3 <- predict(ont_cal1, ensemble=TRUE,
                                    .progress="none",
                                    imp_wts=TRUE, nsim=10,
                                    qvec=NULL))
    expect_equal(dim(pp3),c(14,123,10))
    expect_equal(length(attr(pp3,"imp_wts")),10)
})

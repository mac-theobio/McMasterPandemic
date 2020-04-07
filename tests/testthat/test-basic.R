library(testthat)
library(McMasterPandemic)
library(ggplot2)

context("very basic simulation")

params <- read_params(system.file("params","ICU1.csv",
                                  package="McMasterPandemic"))
state <- make_state(params=params)

test_that("basic examples", {
    expect_is(params,"params_pansim")
    s0 <- run_sim_range(params,state, nt=100)
    expect_is(s0,"data.frame")
    expect_is(state,"state_pansim")
    s1 <- run_sim(params,state,start_date="1-March-2020",end_date="1-Jun-2020")
    expect_is(s1,"pansim")
})

test_that("params methods", {
    expect_equal(summary(params),
                 c(r0 = 0.22783, R0 = 6.51918, Gbar = 12.19064,
                   dbl_time = 3.04243),
                 tolerance=2e-3)
})


time_pars <- data.frame(Date=c("10-Mar-2020","25-Mar-2020"),
                        Symbol=c("beta0","beta0"),
                        Relative_value=c(0.5,0.1))

test_that("time-varying example", {
    resICU_t <- run_sim(params,state,
                        start_date="1-Mar-2020",
                        end_date="1-Jun-2020",
                        params_timevar=time_pars,
                        step_args=list(do_hazard=TRUE))
    expect_is(resICU_t,"pansim")
    plot(resICU_t)
    plot(resICU_t,aggregate=FALSE,log=TRUE,drop_vars=NULL)
})

test_that("time-varying with ndt>1", {
    resICU_t2 <- run_sim(params,state,
                        start_date="1-Mar-2020",
                        end_date="1-June-2020",
                        params_timevar=time_pars,
                        step_args=list(do_hazard=TRUE),
                        ndt=10)
    expect_is(resICU_t2,"pansim")
    plot(resICU_t2)
})

test_that("ndt>1", {
    s2 <- run_sim_range(params,state, nt=100, dt=0.2)
    s3 <- run_sim(params,state, ndt=20,
                  start_date="1-Mar-2020",
                  end_date="20-Mar-2020")
    s3B <- run_sim(params,state,
                  start_date="1-Mar-2020",
                  end_date="20-Mar-2020")
    ## s3C <- run_sim(params,state=unlist(tail(s3B,1)),
    ## start_date="1-Mar-2020",
    ## end_date="20-Mar-2020")
                   
    S3comb <- dplyr::bind_rows(n5=aggregate(s3,pivot=TRUE),
                     n1=aggregate(s3B,pivot=TRUE),
                     .id="model")
    ggplot(S3comb,aes(date,value,colour=var,lty=model))+geom_line() +
        scale_y_log10()
})

test_that("state methods", {
    expect_equal(make_state(N=1,E0=1),
                 structure(c(S = 0, E = 1, Ia = 0, Ip = 0,
                             Im = 0, Is = 0, H = 0, 
                             H2 = 0, ICUs = 0, ICUd = 0,
                             D = 0, R = 0), class = "state_pansim"))
    expect_error(make_state(x=1:5),regexp="must be named")
    expect_warning(make_state(x=c(N=1,E0=1,K=5)),"extra state variables")
})

test_that("calibration", {
    ## test at a tolerance of 0.5%
    s1 <- summary(params)
    s2 <- summary(fix_pars(params,target=c(r=0.23,Gbar=6)))
    expect_equal(s2[["r0"]],0.23,tolerance=0.005)
    expect_equal(s2[["Gbar"]],6,tolerance=0.005)
    s3 <- summary(fix_pars(params,target=c(r=0.23)))
    expect_equal(s3[["r0"]],0.23,tolerance=0.005)
    s4 <- summary(fix_pars(params,target=c(Gbar=6),
                           pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))))
    expect_equal(s4[["Gbar"]],6,tolerance=0.005)
    e1 <- get_evec(params)
    e2 <- get_evec(params,method="analytical")
    expect_equal(e1,e2,tolerance=1e-4)
})

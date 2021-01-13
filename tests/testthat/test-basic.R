library(testthat)
library(McMasterPandemic)
library(ggplot2)

## could run models, create reference models,
##  save("test1","test2","test3","inst/testdata/basic.rda")
##  make target: update_testdata
##  in code: load(system.file("testdata","basic.rda",package="McMasterPandemic"))
context("very basic simulation")

params <- read_params("ICU1.csv")
state <- make_state(params=params,type="ICU1")

test_that("basic examples", {
    expect_is(params,"params_pansim")
    s0 <- run_sim_range(params,state, nt=100)
    expect_is(s0,"data.frame")
    ## original value
    ss <- structure(list(t = 100L, S = 710.7304620948, E = 132.983797059479, 
    Ia = 1280.06036138841, Ip = 9.45207161097982, Im = 2626.01147876609, 
    Is = 47.3684384214429, H = 873.671213383553, H2 = 202.916185137484, 
    ICUs = 626.537390949812, ICUd = 72.4089826500752, D = 3438.87998360902, 
    R = 989978.979634929, foi = 0.00353620556305745),
    state = structure(c(S = 710.7304620948, 
E = 132.983797059479, Ia = 1280.06036138841, Ip = 9.45207161097982, 
Im = 2626.01147876609, Is = 47.3684384214429, H = 873.671213383553, 
H2 = 202.916185137484, ICUs = 626.537390949812, ICUd = 72.4089826500752, 
D = 3438.87998360902, R = 989978.979634929), class = "state_pansim"), row.names = "100", class = "data.frame")
    expect_equal(tail(s0,1), ss,
                 tolerance=1e-8)
    expect_is(state,"state_pansim")
    s1 <- run_sim(params,state,start_date="01-Mar-2020",end_date="01-Jun-2020")
    expect_is(s1,"pansim")
})

test_that("params methods", {
    expect_equal(summary(params),
                 c(r0 = 0.227816539595061, R0 = 6.51800888888889, Gbar = 12.1897401796868, 
                   CFR_gen = 0.0352, dbl_time = 3.0425674175896),
                 tolerance=2e-3)
})


time_pars <- data.frame(Date=c("10-Mar-2020","25-Mar-2020"),
                        Symbol=c("beta0","beta0"),
                        Relative_value=c(0.5,0.1))

test_that("time-varying example", {
    resICU_t <- run_sim(params,state,
                        start_date="01-Mar-2020",
                        end_date="01-Jun-2020",
                        params_timevar=time_pars,
                        step_args=list(do_hazard=TRUE))
    expect_is(resICU_t,"pansim")
    suppressWarnings(print(plot(resICU_t)))
    ## not showing foi because of log / <= 1 filter
    suppressWarnings(print(plot(resICU_t,condense=FALSE,log=TRUE,drop_states=c("t","S","R","E"))))
    ## FIXME: test values!
})

## TO DO: switch dates from anydate()
## test for 'bad date' error
## test multiple param switch


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
    expect_equal(tail(s2,1),
                 tolerance=1e-6,
                 structure(list(t = 100L, S = 999175.599366463, E = 470.938645829701, 
                                Ia = 75.7624206309387, Ip = 26.7668000907004, Im = 128.383139203259, 
                                Is = 5.46913724380533, H = 2.03123004779139, H2 = 0.0400421128254433, 
                                ICUs = 0.371377275540251, ICUd = 0.300078574043836, D = 0.139442448949961, 
                                R = 114.198320079343, foi = 0.000211127356958391), state = structure(c(S = 999175.599366463, 
                                                                                                       E = 470.938645829701, Ia = 75.7624206309387, Ip = 26.7668000907004, 
                                                                                                       Im = 128.383139203259, Is = 5.46913724380533, H = 2.03123004779139, 
                                                                                                       H2 = 0.0400421128254433, ICUs = 0.371377275540251, ICUd = 0.300078574043836, 
                                                                                                       D = 0.139442448949961, R = 114.198320079343), class = "state_pansim"), row.names = "100", class = "data.frame"))

    s3 <- run_sim(params,state, ndt=20,
                  start_date="1-Mar-2020",
                  end_date="20-Mar-2020")
    s3B <- run_sim(params,state,
                  start_date="1-Mar-2020",
                  end_date="20-Mar-2020")
    ## s3C <- run_sim(params,state=unlist(tail(s3B,1)),
    ## start_date="1-Mar-2020",
    ## end_date="20-Mar-2020")
                   
    S3comb <- dplyr::bind_rows(n5=pivot(condense(s3)),
                     n1=pivot(condense(s3B)),
                     .id="model")
    ggplot(S3comb,aes(date,value,colour=var,lty=model))+geom_line() +
        scale_y_log10()
})

test_that("state methods", {
    expect_equal(make_state(N=1,E0=1),
                 structure(c(S = 0, E = 1, Ia = 0, Ip = 0,
                             Im = 0, Is = 0, H = 0,
                             H2 = 0, 
                             ICUs = 0, ICUd = 0,
                             D = 0, R = 0, X=0), class = "state_pansim"))
    expect_error(make_state(x=1:5),regexp="must be named")
    expect_warning(make_state(x=c(N=1,E0=1,K=5)),"extra state variables")
})

test_that("calibration", {
    ## test at a tolerance of 0.5%
    s1 <- summary(params)
    s2 <- summary(p2 <- fix_pars(params,target=c(r=0.23,Gbar=6)))
    expect_equal(s2[["r0"]],0.23,tolerance=0.005)
    expect_equal(s2[["Gbar"]],6,tolerance=0.005)
    ## FIXME: expsim method fails (need to lower uniroot tolerance)
    s3 <- summary(fix_pars(p2,target=c(r=0.4), r_method="rmult"))
    expect_equal(s3[["r0"]],0.4,tolerance=0.005)
    ## adjusting r0 shouldn't change Gbar from previous value
    expect_equal(s3[["Gbar"]], s2[["Gbar"]])
    s4 <- summary(fix_pars(params,target=c(Gbar=6),
                           pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))))
    expect_equal(s4[["Gbar"]],6,tolerance=0.005)
    e1 <- get_evec(params)
    e2 <- get_evec(params,method="analytical")
    expect_equal(e1,e2,tolerance=1e-4)
})

test_that("cum_rep vs rep", {
    params <- update(params,obs_disp=1)
    s0 <- run_sim(params, stoch=c(obs=TRUE,proc=FALSE))
    r <- s0$report
    r[is.na(r)] <- 0
    expect_equal(diff(s0$cumRep), r[-1])
    if (FALSE) {
        ## graphical checking
        matplot(s0[,1],s0[c("report","cumRep")],type="b",log="y")
        op <- par(mfrow=c(1,2))
        plot(report~date,data=s0,type="l")
        plot(cumRep~date,data=s0,type="l")
        par(op)
    }
})

test_that("var-specific obsdisp", {
    params <- update(params,obs_disp=1,obs_disp_I=NA,obs_disp_E=1000)
    set.seed(101)
    s0 <- run_sim(params, stoch=c(obs=TRUE,proc=FALSE))
    plot(s0,keep_states=c("I","E","report"),log=TRUE)
    expect_equal(tail(s0$I,1),16385.5)
    expect_equal(tail(s0$E,1),31791)
})

test_that("mle prediction", {
    load(system.file("testdata","Ontario_basic.rda",package="McMasterPandemic"))
    ## suppress "dropped switch times on final day" warning
    suppressWarnings(test_mle_pred <- predict(Ontario_fit))
    ## hack around test comparison 
    test_mle_pred$var <- unname(test_mle_pred$var)
    ## FIXME: recalculate and save Ontario_basic after conservation fixes
    ## expect_equal(test_mle_pred,mle_prediction)
})


library("McMasterPandemic")
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())

L <- load(system.file("testdata","calib_test.RData",package="McMasterPandemic"))
summary(cparams)
## fill in new parameters for case reports
cparams[["c_prop"]] <- 1/10
cparams[["c_delay_mean"]] <- 5
cparams[["c_delay_cv"]] <- 0.25

cparams[["obs_disp"]] <- 20

## FIXME: thinning interacts with ndt?
set.seed(101)
sim1S <- run_sim(cparams, cstate, start_date="1-Mar-2020",
                 end_date="31-Mar-2020",
                 ndt=10,
                 stoch=c(obs = TRUE, proc = FALSE))

## aggregate/subset simulated data to a short time window (15 Mar - 29 Mar)/
simdat <- (pivot(condense(sim1S))
    %>% filter(date>as.Date("2020-03-15") & date<as.Date("2020-03-29"),
               var %in% c("H","ICU","death","report"))
)

## simdat %>% filter(date==as.Date("2020-03-17"),var=="D")
## sim1S %>% filter(date==as.Date("2020-03-17")) %>% select(D)
## reasonable trajectories
## note H > ICU > D as it should be
plot(sim1S,log=TRUE) + geom_point(data=simdat)

## DRY: condense this some more

## extract H data and set t==0 at beginning
regdatS <- (simdat
    %>% mutate(t0=as.numeric(date-min(date)))
    %>% filter(var=="H")
)
##  fit to H data
g1S <- MASS::glm.nb(value~t0,data=regdatS)

## too slow for now ...
if (FALSE) {
    schoolClose <- "2020-03-17"
    countryClose <- "2020-03-23"
    socialClose <- "2020-03-28"
    bd <- as.Date(c(schoolClose,countryClose,socialClose))

    opt_pars <- list(
        ## these params are part of the main parameter vector: go to run_sim()
        params=c(log_E0=4      ## initial exposed
               , log_beta0=-1  ## initial baseline transmission
                 ## fraction of mild (non-hosp) cases
               , log_mu=log(cparams[["mu"]])
                 ## fraction of incidence reported
                 ## logit_c_prop=qlogis(params[["c_prop"]]),
                 ## fraction of hosp to acute (non-ICU)
               , logit_phi1=qlogis(cparams[["phi1"]])
                 ## fraction of ICU cases dying
                 ## logit_phi2=qlogis(params[["phi2"]])
                 ),
        ## changes in beta at breakpoints
        log_rel_beta0 = rep(-1, length(bd)),
        ## NB dispersion
        log_nb_disp=0)


    ont_all_sub <- (ont_all
        %>% mutate_at("var",trans_state_vars)
        %>% filter(var %in% c("H","ICU","death","report"))
    )

    params <- fix_pars(read_params("ICU1.csv")
                     , target=c(Gbar=6)
                     , pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
                       )
    params[["N"]] <- 14.57e6  ## reset pop to Ontario

    system.time(cc1 <- calibrate(data=ont_all_sub
                               , base_params=params
                               , opt_pars = opt_pars
                               , time_args=list(break_dates = bd)
                                 )
                )

    system.time(cc2 <- calibrate(data=ont_all_sub
                               , base_params=params
                               , opt_pars = opt_pars
                               , time_args=list(break_dates = bd)
                               , sim_args=list(use_ode=TRUE)
                                 )
                )

    repdata <- ont_all_sub %>% filter(var=="report")
    cc3 <- calibrate_comb(data=repdata
                        , params=params
                        , use_phenomhet=TRUE
                        , use_spline=FALSE
                          )
} ## end slow stuff

#### make sure single-variable fit works OK

p1 <- fix_pars(read_params("ICU1.csv"))
p2 <- update(p1, obs_disp=1, proc_disp=0, zeta=5)
set.seed(101)
r1 <- run_sim(p2, stoch=c(obs=TRUE, proc=TRUE), end_date="2020-05-31")

dd_r <- r1 %>% select(date,report) %>% pivot() %>% na.omit()

p3 <- update(p2, E0=exp(4), beta0=exp(-1))  ## set parameters to *original* starting values from get_opt_pars
c_r2 <- calibrate_comb(params=p3,
                       use_phenomhet=FALSE,
                       debug_plot=FALSE,
                       data=dd_r, use_DEoptim=FALSE,
                       use_spline=FALSE)

plot(c_r2, data=dd_r)
## list(params = c(E0 = 0.969447127312371,
##                 beta0 = 0.999559822048325
##                 ),
##      nb_disp = c(report = 1.11651207216957),
##      time_beta = numeric(0)),

## CHANGED again with conservation fix
## ref_val <- list(params = c(E0 = 4.35599429735943, beta0 = 0.804378506840628),
##                 nb_disp = c(report = 1.13530065646639), time_beta = numeric(0))

## CHANGED but these look better ... ?
## ref_val <- list(params = c(E0 = 2.22166438860786, beta0 = 0.873467646391076),
##                 nb_disp = c(report = 0.996186838808113), time_beta = numeric(0))

## CHANGED again (X/hosp accumulator)
ref_val <- list(params = c(E0 = 63.3150461392819, beta0 = 0.649997557506806),
                nb_disp = c(report = 0.495604121823216), time_beta = numeric(0))

## ref_val <- list(params = c(E0 = 2.22166438860786, beta0 = 0.873467646391076),
##                 nb_disp = c(report = 0.996186838808113), time_beta = numeric(0))

## broken in 54095c?
print(ref_val)
print(coef(c_r2, "fitted"))
stopifnot(all.equal(coef(c_r2,"fitted"),
                    ref_val, 
                    tolerance=1e-6))

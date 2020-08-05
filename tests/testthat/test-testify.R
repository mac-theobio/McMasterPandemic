library(McMasterPandemic)
library(testthat)

context("testify")

pp <- read_params("PHAC_testify.csv")

## Making states and expanding states
state <- make_state(params=pp)
state_testified <- expand_stateval(state)

lfun <- function(x,newstates=4,fixedstates=3) {
	newstates*(length(x)-fixedstates)+fixedstates
}

## fixedstates are D,R,X
lfun(state)

test_that("testified states make sense", {
    expect_equal(length(state_testified), lfun(state)+2)
    expect_equal(sort(unique(gsub("_.*$","",names(state_testified)))),
                 c(sort(c("N","P",names(state)))))
})

## Making beta_vec wtr states (infectious compartments only)
beta_vec0 <- make_betavec(state,pp,full=FALSE)
beta_vec0_testified <- make_betavec(state_testified,pp,full=FALSE)
## full
beta_vec <- make_betavec(state,pp)
beta_vec_testified <- make_betavec(state_testified,pp)

test_that("testified betas make sense", {
    ## test names of testified beta against I_
    expect_equal(names(beta_vec0_testified), grep("^I[[:lower:]]",names(state_testified),value=TRUE))
    expect_equal(unname(beta_vec0_testified[grepl("_u$",names(beta_vec0_testified))]),
                 unname(beta_vec0))
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
g
## Updating FOI
test_that("FOI doesn't change", {
    pp2 <- update(pp,iso_p=0)
    expect_equal(update_foi(state,pp2,beta_vec),
                 update_foi(state_testified,pp2,beta_vec_testified))
})

test_that("condensation is OK", {
    sim0_testified_condensed <- run_sim(params = pp,
                                        ratemat_args = list(testify=TRUE))
    expect_equal(names(sim0_testified_condensed),
                 c("date", "S", "E", "I", "H", "hosp", "ICU", "R",
                   "death", "D", 
                   "negtest", "N", "postest", "P",
                   "foi", "incidence", "report", 
                   "cumRep"))
})


test_that("time-varying test intensity", {
    pt <- data.frame(Date=as.Date(c("2020-04-01","2020-04-15")),
                     Symbol=rep("testing_intensity",2),
                     Relative_value=c(0,4))
    pp[["testing_intensity"]] <- 0.002
    sim0_testified_timevar <- run_sim(params = pp,
                                      ratemat_args = list(testify=TRUE),
                                      params_timevar=pt)
    ## library(directlabels)
    p1 <- plot(sim0_testified_timevar,log=TRUE,log_lwr=1e-7)
    ## direct.label(p1,method="last.bumpup")
    pvars <- c("N","P","negtest","postest")
    expect_equal(unlist(tail(sim0_testified_timevar[,pvars],1)),
                 c(N = 35338.9675996009, P = 10.2554061586176,
                   negtest = 775.310922278855, 
                   postest = 3.97882806061493))
    ##
    pp_noinf <- update(pp,beta0=0,E0=0)
    ## inf etc. is non-zero because we start with E0>0 (but trivial)
    sim0_noinf <- run_sim(params = pp_noinf,
                          ratemat_args = list(testify=TRUE))
    sim0_noinf_all <- run_sim(params = pp_noinf,
                              ratemat_args = list(testify=TRUE),
                              condense_args = list(keep_all=TRUE,
                                               add_reports=FALSE))
    ## negtest *should* converge on
    pp[["testing_intensity"]]*pp[["N"]]
    tail(sim0_noinf,1)[c("N","negtest")]
    head(sim0_noinf)
    head(sim0_noinf_all)
    plot(sim0_noinf,log=TRUE)
    nm <- c("S","E","I","H","ICU")
    sum(tail(sim0_testified_timevar[,nm],1))*pp[["testing_intensity"]]
    ## testing intensity = 0.002
    ## omega = 0.1 (wating time=10 days)
    ## du/dt = -i*u+omega*n
    ## dn/dt =  i*u-omega*n
    ## steady state of {u,n}? n = (i/omega)*u ~ (i/omega)*S
    ## 
    ## 
    ## dN/dt = omega*n -> ~ i*S
    ##
    
    sim0_noinf_all %>% dplyr::select(S_u,S_n,N) %>%
        dplyr::mutate(negtest=c(NA,diff(N))) %>% tail()
    
    plot(sim0_noinf_all, log=TRUE)
})

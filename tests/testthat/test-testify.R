library(McMasterPandemic)
library(testthat)

## squash call so we ignore diffs when using all.equal()
uncall <- function(x) {
    attr(x,"call") <- NULL
    return(x)
}

pp <- read_params("PHAC_testify.csv")
## Making states and expanding states
state <- make_state(params=pp,testify=FALSE)
state_testified <- expand_stateval_testing(state, params=pp, method="untested")

## global variables
## print(non_expanded_states)
## print(test_extensions)
## print(test_accumulators)

lfun <- function(x,n_exts=length(test_extensions),n_fixed=length(non_expanded_states)) {
    n_exts*(length(x)-n_fixed)+n_fixed
}

lfun(state)

test_that("testified states make sense", {
    expect_equal(length(state_testified), lfun(state)+2)
    expect_equal(sort(unique(gsub("_.*$","",names(state_testified)))),
                 c(sort(c("N","P",names(state)))))
})

test_that("make_state from scratch", {
    expect_equal(names(make_state(params=pp)),
                 names(state_testified))
})

test_that("make_ratemat from scratch (ignore testify in state)", {
    s0 <- make_state(params=pp, testify=FALSE)
    M0 <- make_ratemat(s0, pp)
    s1 <- make_state(params=pp)
    M1 <- make_ratemat(state=s1,params=pp)
    ## FIXME: why does testify/untestify swap order of D and R?
    ## should be harmless but ...
    M1 <- M1[rownames(M0),colnames(M0)]
    expect_equal(dim(M0),dim(M1))
})

## Making beta_vec wtr states (infectious compartments only)
beta_vec0 <- make_beta(state,pp,full=FALSE)
beta_vec0_testified <- make_beta(state_testified,pp,full=FALSE)
## full
beta_vec <- make_beta(state,pp)
beta_vec_testified <- make_beta(state_testified,pp)

test_that("testified betas make sense", {
    ## test names of testified beta against I_
    expect_equal(names(beta_vec0_testified), grep("^I[[:lower:]]",names(state_testified),value=TRUE))
    expect_equal(unname(beta_vec0_testified[grepl("_u$",names(beta_vec0_testified))]),
                 unname(beta_vec0))
})

test_that("pos test vec check", {
    vn <- setdiff(names(state),non_expanded_states)
    badp <- c(Pxx=1)
    expect_error(make_test_posvec(badp,vn), "vector names should match")
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

## Updating FOI
test_that("FOI doesn't change (for 'untested' expansion)", {
    pp2 <- update(pp,iso_t=0)
    expect_equal(update_foi(state,pp2,beta_vec),
                 update_foi(state_testified,pp2,beta_vec_testified))
})

sim0_testified_condensed <- run_sim(params = pp,
                                    ## specify testing time to avoid warning
                                    ratemat_args = list(testing_time="sample"))

make_state(pp[["N"]], pp[["E0"]], params=pp)
test_that("obsolete testify spec", {
    expect_warning(run_sim(params = pp,
                           ratemat_args = list(testify=TRUE)),
                   "no longer needs to be passed")
})

test_that("testing time default", {
          expect_warning(sim0_tc2 <-  run_sim(params = pp),
                         "setting testing time to 'sample'")
          expect_equal(uncall(sim0_tc2), uncall(sim0_testified_condensed))
})

test_that("condensation is OK", {
    expect_equal(names(sim0_testified_condensed),
                 c("date", "S", "E", "I", "H", "ICU", "R",
                   "hosp", "X", "death", "D",
                   "negtest", "N", "postest", "P",
                   "foi", "incidence", "report",
                   "cumRep"))
 })

test_that("time-varying test intensity", {
     pt <- data.frame(Date=as.Date(c("2020-04-01","2020-04-15")),
                      Symbol=rep("testing_intensity",2),
                      Relative_value=c(0.1,4))
     ## don't reduce testing to 0 - it will break things!
     pp[["testing_intensity"]] <- 0.002
     sim0_testified_timevar <- run_sim(params = pp,
                                       ratemat_args = list(testing_time="sample"),
                                       params_timevar=pt)
#     ## library(directlabels)
     p1 <- plot(sim0_testified_timevar,log=TRUE,log_lwr=1e-7)
#     ## direct.label(p1,method="last.bumpup")
     pvars <- c("N","P","negtest","postest")
#     expect_equal(unlist(tail(sim0_testified_timevar[,pvars],1)),
#                  c(N = 75642.0699087884, P = 55.0789129220698,
#                    negtest = 4971.36249153456,
#                    postest = 21.6243455422309))
})
#
test_that("testing with susceptibles only", {
     pp[["testing_intensity"]] <- 0.002
     pp_noinf <- update(pp,beta0=0,E0=0)  ## no transmission, no infected people
     ## suppress warning about "initial values too small for rounding"
     sim0_noinf <- suppressWarnings(run_sim(params = pp_noinf,
                           ratemat_args = list(testing_time="sample"),
                           end_date="2021-01-01"))
     ## negtest *should* converge on:
     expected_negtest <- with(as.list(pp),omega/(omega+testing_intensity)*testing_intensity*N)
     ## expect_equal(tail(sim0_noinf[["negtest"]],1),expected_negtest)
})

plotfun <- function(L) {
     if (require(dplyr) && require(ggplot2) && require(purrr) && require(tidyr)) {
         dd <- (map_dfr(L,
                        ~(select(.,c(date,negtest,postest,report))),.id="type")
             %>% pivot_longer(cols=c(negtest,postest,report))
         )
         return(suppressWarnings(ggplot(dd,
                                        aes(date,value,colour=name,linetype=type))
                                 + geom_line() + scale_y_log10())
                )
     }
 }

test_that("testify + sampling time", {
     rvars <- c("N","P")
     sim0_testified_report <- run_sim(params = pp,
                                      ratemat_args = list(testing_time="report"))
     res_report <- tail(sim0_testified_report[rvars],1)
     sim0_testified_sample <- run_sim(params = pp,
                                      ratemat_args = list(testing_time="sample"))
     res_sample <- tail(sim0_testified_sample[rvars],1)
     expect(!(all(res_report==res_sample)), "results don't differ according to testing time")
#     expect_equal(unlist(res_sample),c(N = 78578.6760706341, P = 107.200806051815))
     gg1 <- plotfun(list(report=sim0_testified_report,sample=sim0_testified_sample))
})
#
 test_that("testify + alternative weights models", {
     rvars <- c("N","P")
     ppw0 <- pp[!grepl("^W",names(pp))]  ## remove all of the regular W-parameters
     class(ppw0) <- "params_pansim"  ## restore class (lost after indexing) so we can use update()
     ppw1 <- update(ppw0,W_asymp=0.2)
     ppw2 <- update(ppw0,c(W_asymp=0.2, W_severe=1.5))
     sim0_w0 <- run_sim(params = pp,
                        ratemat_args = list(testing_time="report"))
     sim0_w1 <- run_sim(params = ppw1,
                        ratemat_args = list(testing_time="report"))
     sim0_w2 <- run_sim(params = ppw2,
                        ratemat_args = list(testing_time="report"))
     gg1 <- plotfun(list(w0=sim0_w0,w1=sim0_w1,w2=sim0_w2))
     ref_val <- structure(list(N = c(62097.0986898706, 61789.5190257482, 61781.6457354169),
                               P = c(13.9744057111373, 16.3878085964286, 16.3569118395614)),
                          row.names = c("43", "431", "432"), class = c("pansim", "data.frame"))
     ## these are *almost* identical (is that what I expected?)
#     expect_equal(ref_val,
#         rbind(tail(sim0_w0[rvars],1),
#           tail(sim0_w1[rvars],1),
#           tail(sim0_w2[rvars],1)))
 })


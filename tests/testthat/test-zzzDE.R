library(McMasterPandemic)
library(testthat)
library(parallel)

test_that("onLoad options", {
    ## expect no errors while loading
    expect(!is.null(options()$MP_badsum_action), "MP_badsum_action not set!")
    expect(!is.null(options()$MP_badsum_tol), "MP_badsum_tol not set!")
    ## expect error for random variable
    expect(
        is.null(options()$random_var_for_testing_please_do_not_name_a_variable_this),
        "This should never happen, random_var_for_testing_please_do_not_name_a_variable_this is not NULL"
    )
})

context("DEoptim")

test_level <- if (nzchar(s <- Sys.getenv("MACPAN_TEST_LEVEL")) &&
    is.finite(s <- as.numeric(s))) {
    s
} else {
    1
}
test_level = 2

library(dplyr)
params <- fix_pars(read_params("ICU1.csv"))
opt_pars <- list(
    params = c(
        log_E0 = 4, log_beta0 = -1,
        logit_mu = plogis(params[["mu"]]),
        logit_phi1 = qlogis(params[["phi1"]])
    ),
    logit_rel_beta0 = c(-1, -1),
    log_nb_disp = NULL
)
dd <- (ont_all %>% trans_state_vars() %>% filter(var %in% c("report", "death", "H")))

## skip on Travis/GitHub Actions
if (Sys.getenv("TRAVIS") != "true" ## not on Travis
&&
    Sys.getenv("CI_WORKFLOW") == "" ## not on GH Actions (use NOT_CRAN?)
&&
    test_level > 1) {
    suppressWarnings(cal1_DE <- calibrate(
        data = dd,
        base_params = params,
        opt_pars = opt_pars,
        use_DEoptim = TRUE,
        DE_cores = 1,
        DE_args = list(control = list(itermax = 5, trace = FALSE)),
        time_args = list(params_timevar = NULL)
    ))

    de <- attr(cal1_DE, "de")

    test_that("DE components are present", {
        expect_is(attr(cal1_DE, "de_time"), "proc_time")
        expect_is(de, "DEoptim")
        expect(!is.null(de$member$Sigma), failure_message = "Sigma component missing")
    })


    suppressWarnings(cal2_DE <- calibrate_comb(
        data = dd,
        params = params,
        use_DEoptim = TRUE,
        DE_cores = 1,
        DE_args = list(control = list(itermax = 5, trace = FALSE))
    ))
} ## skip slow tests

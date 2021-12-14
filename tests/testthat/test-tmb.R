Sys.setenv(R_TESTS="")

library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)

## TODO/FIXME/BACKGROUND --------------------------------------------
## rename/repurpose this once the tmb engine is properly added
## to the package. for now refer to a particular version
## of the spec and store the C++ files in the
## appropriate directory within inst. also for now, assume that this
## test file is being run inside the directory in which it should
## live within the package structure
## (i.e. inside McMasterPandemic/tests/testthat). this should
## also work for `make pkgtest` because tests are run in the
## directory of the test file (as far as I can tell), but it could
## be fragile in other contexts? in any case, once we have
## the tmb-refactored code in the package, we will use the
## installed dlls and R-side utilities but retain the structure of
## these tests by spec version. possibly we can roll-up all up all
## spec versions with the same major version together (i.e. those
## that are backwards compatible).
## see https://canmod.net/misc/flex_specs for more on
## spec versioning.
## ------------------------------------------------------------------

test_that("spec v0.0.1 rate matrices match make_ratemat", {
    set_spec_version("0.0.1", "../../inst/tmb/")

    params <- read_params("ICU1.csv")
    state <- McMasterPandemic::make_state(params = params)
    M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)
    test_model <- (
        init_model(params, state)
            %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
            %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
            %>% add_rate("Ia", "R", ~ (gamma_a))
            %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
            %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
            %>% add_rate("Im", "R", ~ (gamma_m))
            %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
            %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
            %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
            %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
            %>% add_rate("ICUs", "H2", ~ (psi1))
            %>% add_rate("ICUd", "D", ~ (psi2))
            %>% add_rate("H2", "R", ~ (psi3))
            %>% add_rate("H", "R", ~ (rho))
            %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
            %>% add_rate("S", "E", ~
            (Ia) * (beta0) * (1 / N) * (Ca) +
                (Ip) * (beta0) * (1 / N) * (Cp) +
                (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
                (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
            %>% add_tmb_indices()
    )

    dd <- tmb_fun(test_model)

    tmb_sparse_ratemat <- dd$report()$ratemat
    r_dense_ratemat <- as.matrix(M)
    tmb_dense_ratemat <- as.matrix(tmb_sparse_ratemat)
    r_sparse_ratemat <- M
    dimnames(tmb_sparse_ratemat) <- dimnames(r_sparse_ratemat) <- dimnames(tmb_dense_ratemat) <- dimnames(r_dense_ratemat) <- NULL

    expect_equal(
        class(tmb_sparse_ratemat),
        class(r_sparse_ratemat)
    )

    expect_equal(
        tmb_sparse_ratemat@Dim,
        r_sparse_ratemat@Dim
    )

    expect_equal(
        tmb_dense_ratemat,
        r_dense_ratemat
    )
})

test_that("spec v0.0.2 simulations match run_sim", {

    set_spec_version("0.0.2", "../../inst/tmb/")
    params <- read_params("ICU1.csv")
    state <- McMasterPandemic::make_state(params = params)
    test_model <- (
        init_model(params, state = state)
            %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
            %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
            %>% add_rate("Ia", "R", ~ (gamma_a))
            %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
            %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
            %>% add_rate("Im", "R", ~ (gamma_m))
            %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
            %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
            %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
            %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
            %>% add_rate("ICUs", "H2", ~ (psi1))
            %>% add_rate("ICUd", "D", ~ (psi2))
            %>% add_rate("H2", "R", ~ (psi3))
            %>% add_rate("H", "R", ~ (rho))
            %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
            %>% add_rate("S", "E", ~
            (Ia) * (beta0) * (1 / N) * (Ca) +
                (Ip) * (beta0) * (1 / N) * (Cp) +
                (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
                (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
            %>% add_parallel_accumulators(c("X", "N", "P", "V"))
            %>% add_tmb_indices()
    )

    dd <- tmb_fun(test_model)

    tmb_traj <- (test_model$state
        %>% c(dd$report()$concatenated_state_vector)
        %>% matrix(length(state), 4, dimnames = list(names(state), 1:4))
        %>% t()
        %>% as.data.frame()
    )
    r_traj <- run_sim_range(
        params = params, state = test_model$state, nt = 4,
        step_args = list(do_hazard = FALSE)
    )[, names(test_model$state)]

    expect_equal(tmb_traj, r_traj)
})

test_that("spec v0.0.4 simulations with time varying parameters match run_sim", {
    set_spec_version("0.0.4", "../../inst/tmb/")
    params <- read_params("ICU1.csv")
    tv_dat <- data.frame(
        Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
        Symbol = c("beta0", "beta0", "beta0"),
        Value = c(0.5, 0.1, 0.05),
        Type = c("rel_prev", "rel_orig", "rel_prev")
    )

    test_model <- (init_model(
        params,
        state = make_state(params = params),
        start_date = "2021-09-09", end_date = "2021-10-09",
        params_timevar = tv_dat
    )
        %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
        %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
        %>% add_rate("Ia", "R", ~ (gamma_a))
        %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
        %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
        %>% add_rate("Im", "R", ~ (gamma_m))
        %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
        %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
        %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
        %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
        %>% add_rate("ICUs", "H2", ~ (psi1))
        %>% add_rate("ICUd", "D", ~ (psi2))
        %>% add_rate("H2", "R", ~ (psi3))
        %>% add_rate("H", "R", ~ (rho))
        %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
        %>% add_rate("S", "E", ~
        (Ia) * (beta0) * (1 / N) * (Ca) +
            (Ip) * (beta0) * (1 / N) * (Cp) +
            (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
            (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
        %>% add_parallel_accumulators(c("X", "N", "P", "V"))
        %>% add_tmb_indices()
    )

    dd <- tmb_fun(test_model)

    tmb_traj <- (test_model$state
        %>% c(dd$report()$concatenated_state_vector)
        %>% matrix(length(test_model$state), test_model$iters + 1,
            dimnames = list(names(test_model$state), 1:(test_model$iters + 1))
        )
        %>% t()
        %>% as.data.frame()
    )

    r_traj <- run_sim(
        params = params,
        state = test_model$state,
        start_date = test_model$start_date,
        end_date = test_model$end_date,
        params_timevar = tv_dat,
        condense = FALSE,
        step_args = list(do_hazard = FALSE)
    )[, names(test_model$state)] %>%
        as.data.frame()

    row.names(r_traj) <- row.names(tmb_traj) <- NULL

    expect_equal(tmb_traj, r_traj)
})

test_that("spec v0.0.5 simulations with hazard steps match run_sim, and autodiff is working", {

    set_spec_version("0.0.5", "../../inst/tmb/")
    params <- read_params("ICU1.csv")
    state = make_state(params = params)

    tv_dat <- data.frame(
        Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
        Symbol = c("beta0", "beta0", "beta0"),
        Value = c(0.5, 0.1, 0.05),
        Type = c("rel_prev", "rel_orig", "rel_prev")
    )


    test_model <- (init_model(
        params,
        state = make_state(params = params),
        start_date = "2021-09-09", end_date = "2021-10-09",
        params_timevar = tv_dat
    )
        %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
        %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
        %>% add_rate("Ia", "R", ~ (gamma_a))
        %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
        %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
        %>% add_rate("Im", "R", ~ (gamma_m))
        %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
        %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
        %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
        %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
        %>% add_rate("ICUs", "H2", ~ (psi1))
        %>% add_rate("ICUd", "D", ~ (psi2))
        %>% add_rate("H2", "R", ~ (psi3))
        %>% add_rate("H", "R", ~ (rho))
        %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
        %>% add_rate("S", "E", ~
        (Ia) * (beta0) * (1 / N) * (Ca) +
            (Ip) * (beta0) * (1 / N) * (Cp) +
            (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
            (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
        %>% add_parallel_accumulators(c("X", "N", "P", "V"))
        %>% add_tmb_indices()
    )

    dd <- tmb_fun(test_model)

    tmb_traj <- (dd$report()$concatenated_state_vector
        %>% matrix(length(test_model$state), test_model$iters + 1,
            dimnames = list(names(test_model$state), 1:(test_model$iters + 1))
        )
        %>% t()
        %>% as.data.frame()
    )

    r_traj <- run_sim(
        params = params,
        state = test_model$state,
        start_date = test_model$start_date,
        end_date = test_model$end_date,
        params_timevar = tv_dat,
        condense = FALSE,
        step_args = list(do_hazard = TRUE)
    )[, names(test_model$state)] %>%
        as.data.frame()

    row.names(r_traj) <- row.names(tmb_traj) <- NULL

    # simulated state trajectories are equal
    expect_equal(tmb_traj, r_traj)

    # mock objective function is correct
    expect_equal(sum(tmb_traj[31, ]), dd$fn(dd$par))

    # tmb-computed gradient equals numerical gradient

    ## numerical differentiation settings:
    ## used defaults from `?grad` with one exception: d = 0.1, not 0.0001.
    ## i don't understand why this helps, but it allows us to take the
    ## tolerance from 1e-4 down to 1e-6
    numerical_deriv_args <-
        list(
            eps = 1e-4, d = 0.1,
            zero.tol = sqrt(.Machine$double.eps / 7e-7), r = 4, v = 2,
            show.details = FALSE
        )


    numeric_gradient <- numDeriv::grad(dd$fn, dd$par,
        method.args = numerical_deriv_args
    )
    tmb_gradient <- dd$gr(dd$par)
    attributes(numeric_gradient) <- attributes(tmb_gradient) <- NULL
    expect_equal(tmb_gradient, numeric_gradient,
        tolerance = 1e-5
    )


    # use_flex flag does not change results
    tmb_sim <- run_sim(
        params, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        params_timevar = tv_dat,
        condense = FALSE,
        step_args = list(do_hazard = TRUE),
        use_flex = TRUE
    )
    r_sim <- run_sim(
        params = params, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        params_timevar = tv_dat,
        condense = FALSE,
        step_args = list(do_hazard = TRUE),
        use_flex = FALSE
    )
    compare_sims(r_sim, tmb_sim)
})

test_that('spec v0.0.6 time-varying parameters are correctly updated on C++ side', {
    set_spec_version("0.0.6", "../../inst/tmb/")
    params <- read_params("ICU1.csv")
    state <- make_state(params = params)
    M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

    tv_dat <- data.frame(
        Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
        Symbol = c("beta0", "beta0", "beta0"),
        Value = c(0.5, 0.1, 0.05),
        Type = c("rel_prev", "rel_orig", "rel_prev")
    )

    mm = make_base_model(
        params, state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        params_timevar = tv_dat,
        step_args = list(do_hazard = TRUE)
    )

    # change beta0, which is time-varying, so that
    # we can check that the parameter updates are
    # happening correctly in the C++ side
    test_pars = params
    test_pars[1] = 3

    tmb_sim <- run_sim(
        params = test_pars, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        params_timevar = tv_dat,
        condense = FALSE,
        step_args = list(do_hazard = TRUE, flexmodel = mm),
        use_flex = TRUE
    )

    r_sim <- run_sim(
        params = test_pars, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        params_timevar = tv_dat,
        condense = FALSE,
        step_args = list(do_hazard = TRUE),
        use_flex = FALSE
    )
    compare_sims(r_sim, tmb_sim)
})

test_that('spec v0.0.6 that it remains ok to _not_ use time-varying parameters', {
    set_spec_version("0.0.6", "../../inst/tmb/")
    params <- read_params("ICU1.csv")
    state <- make_state(params = params)
    M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

    tmb_sim <- run_sim(
        params = params, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        condense = TRUE,
        step_args = list(do_hazard = TRUE),
        use_flex = TRUE
    )

    r_sim <- run_sim(
        params = params, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        condense = TRUE,
        step_args = list(do_hazard = TRUE),
        use_flex = FALSE
    )
    compare_sims(r_sim, tmb_sim)
})

test_that('spec v0.1.1 tmb outflow can be set to match exponential simulation', {
    reset_spec_version()

    params <- read_params("ICU1.csv")

    # modify parameters and state for eigenvector calculation ('by hand')
    state <- make_state(params = params)[1:12]
    params[['N']] = 1
    params[['E0']] = 1e-5
    state[] = 0
    state[['S']] = 1 - params[['E0']]
    state[['E']] =     params[['E0']]
    iters = 100

    r_sim = run_sim_range(params, state, nt = iters,
                          step_args = list(do_hazard = FALSE,
                                           do_exponential = TRUE))

    start_date = ymd(20000101)
    model <- (init_model(params, state,
                         start_date = start_date,
                         end_date = start_date + days(iters - 1),
                         do_hazard = FALSE,
                         do_make_state = FALSE)
              %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
              %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
              %>% add_rate("Ia", "R", ~ (gamma_a))
              %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
              %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
              %>% add_rate("Im", "R", ~ (gamma_m))
              %>% add_rate("Is", "H", ~
                               (1 - nonhosp_mort) * (phi1) * (gamma_s))
              %>% add_rate("Is", "ICUs", ~
                               (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
              %>% add_rate("Is", "ICUd", ~
                               (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
              %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
              %>% add_rate("ICUs", "H2", ~ (psi1))
              %>% add_rate("ICUd", "D", ~ (psi2))
              %>% add_rate("H2", "R", ~ (psi3))
              %>% add_rate("H", "R", ~ (rho))
              %>% add_rate("S", "E", ~
                               (Ia) * (beta0) * (1 / N) * (Ca) +
                               (Ip) * (beta0) * (1 / N) * (Cp) +
                               (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
                               (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
              %>% add_outflow("^S$", "^S$")
              %>% add_outflow("^(E|I|H|ICU|D|R)", "^(S|E|I|H|ICU|D|R)")
              %>% add_tmb_indices
    )

    # trim off S, D, & R -- for eigenvector calculation 'by hand'
    r_final_state = unlist(r_sim[100, names(state)])[2:10]
    tmb_final_state = c(final_state_vector(model)[2:10])

    expect_equal(norm_vec(r_final_state), norm_vec(tmb_final_state))
})

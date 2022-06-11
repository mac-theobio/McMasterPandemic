Sys.setenv(R_TESTS="")

library(testthat)
library(McMasterPandemic)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)
library(here)

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
pkg_home = here()
inst_tmb = file.path(pkg_home, 'inst/tmb')

test_that("spec v0.0.1 rate matrices match make_ratemat", {
    set_spec_version("0.0.1", system.file('tmb', package = 'McMasterPandemic'))
    options(MP_force_dgTMatrix = TRUE)
    params <- read_params("ICU1.csv")
    state <- McMasterPandemic::make_state(params = params)
    M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)
    test_model <- (
        flexmodel(params, state)
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
            %>% update_tmb_indices()
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

    set_spec_version("0.0.2", system.file('tmb', package = 'McMasterPandemic'))
    params <- read_params("ICU1.csv")
    state <- McMasterPandemic::make_state(params = params)
    test_model <- (
        flexmodel(params, state = state)
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
            %>% update_tmb_indices()
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
    set_spec_version("0.0.4", system.file('tmb', package = 'McMasterPandemic'))
    params <- read_params("ICU1.csv")
    tv_dat <- data.frame(
        Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
        Symbol = c("beta0", "beta0", "beta0"),
        Value = c(0.5, 0.1, 0.05),
        Type = c("rel_prev", "rel_orig", "rel_prev")
    )

    test_model <- (flexmodel(
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
        %>% update_tmb_indices()
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

    set_spec_version("0.0.5", system.file('tmb', package = 'McMasterPandemic'))
    r_tmb_comparable()
    params <- read_params("ICU1.csv")
    state = make_state(params = params)

    tv_dat <- data.frame(
        Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
        Symbol = c("beta0", "beta0", "beta0"),
        Value = c(0.5, 0.1, 0.05),
        Type = c("rel_prev", "rel_orig", "rel_prev")
    )


    test_model <- (flexmodel(
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
        %>% update_tmb_indices()
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
        use_flex = TRUE,
        flexmodel = make_base_model(
          params = params,
          state = state,
          start_date = "2021-09-10",
          end_date = "2021-10-10",
          params_timevar = tv_dat,
          do_make_state = FALSE,
          do_hazard = TRUE
        )
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
    set_spec_version("0.0.6", system.file('tmb', package = 'McMasterPandemic'))
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
        use_flex = TRUE,
        flexmodel = mm
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
    set_spec_version("0.0.6", system.file('tmb', package = 'McMasterPandemic'))
    params <- read_params("ICU1.csv")
    state <- make_state(params = params)
    M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

    tmb_sim <- run_sim(
        params = params, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        condense = TRUE,
        step_args = list(do_hazard = TRUE),
        use_flex = TRUE,
        flexmodel = make_base_model(
          params = params,
          state = state,
          start_date = "2021-09-10",
          end_date = "2021-10-10",
          do_make_state = FALSE,
          do_hazard = TRUE
        )
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

    params <- read_params("PHAC.csv")

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
    model <- (flexmodel(params, state,
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
              %>% update_tmb_indices
    )

    # trim off S, D, & R -- for eigenvector calculation 'by hand'
    r_final_state = unlist(r_sim[100, names(state)])[2:10]
    tmb_final_state = c(final_state_vector(model)[2:10])
    norm_vec = McMasterPandemic:::norm_vec
    expect_equal(norm_vec(r_final_state), norm_vec(tmb_final_state))
})

test_that("flex models made with null state can be used", {
    reset_spec_version()
    r_tmb_comparable()
    params = read_params("PHAC.csv")
    model = make_base_model(
        params = params,
        start_date = "1900-01-01",
        end_date = "1900-02-01")
    tmb_sims = run_sim(
        params = params,
        start_date = model$start_date,
        end_date = model$end_date,
        flexmodel = model
    )
    r_sims = run_sim(
        params = params,
        start_date = model$start_date,
        end_date = model$end_date
    )
    compare_sims(r_sims, tmb_sims, na_is_zero = TRUE)
})

test_that("simple sir models produce correct simulations", {
    S0 = 20000
    I0 = 100
    R0 = 0
    sir_model = (
        flexmodel(
            params = c(
                N = S0 + I0,
                gamma = 0.06, # per-infected recovery rate
                beta = 0.15   # per-contact transmission rate
            ),
            state = c(S = S0, I = I0, R = R0),
            start_date = "2000-01-01",
            end_date = "2000-05-01",
            do_hazard = FALSE,
            do_make_state = FALSE
        )
        %>% add_rate("S", "I", ~ (1/N) * (beta) * (I))
        %>% add_rate("I", "R", ~ (gamma))
        %>% add_outflow(".+", ".+")
        %>% update_tmb_indices
    )

    i <- 1
    lenSim <- 122
    S <-  E <-  I <-  R <- D <- numeric(lenSim)
    S[1] <- 20000
    I[1] <- 100
    N <- S[1] + I[1]
    R[1] <- 0
    gamma <- 0.06
    beta<- 0.15

    while (i <= lenSim - 1){
        S[i+1] <- S[i] - (beta) * (S[i]) * (I[i])/N
        I[i+1] <- I[i] - gamma * I[i] + (beta) * (S[i]) * (I[i])/N
        R[i+1] <- R[i]  + gamma * I[i]
        i <- i + 1
    }

    sim <- data.frame("S" = S, "I" = I, "R" = R, "t" = 1:lenSim)
    macpan_sim <- simulate_state_vector(sir_model)

    expect_equal(macpan_sim$S, sim$S)
    expect_equal(macpan_sim$I, sim$I)
    expect_equal(macpan_sim$R, sim$R)
})

test_that("one may specify different rates for the same flow", {
    reset_spec_version()
    tmb_mode()
    options(MP_warn_repeated_rates = TRUE)
    model_rep = model_one = (flexmodel(
            params = c(alpha = 0.1),
            state = c(A = 100, B = 0),
            start_date = "2000-01-01",
            end_date = "2000-02-01",
            do_hazard = FALSE, do_make_state = FALSE
        )
          %>% add_rate("A", "B", ~ (alpha))
          %>% add_outflow()
    )
    expect_warning(model_rep <- (model_rep
      %>% add_rate("A", "B", ~ (alpha))
      %>% update_tmb_indices()
    ))
    model_one$params = c(alpha = 0.2)
    model_one = update_tmb_indices(model_one)
    expect_equal(
        simulate_state_vector(model_one),
        simulate_state_vector(model_rep))
})

test_that("an informative error is returned if variables are missing", {
    reset_spec_version()
    tmb_mode()
    tv = data.frame(
        Date = "2000-02-01",
        Symbol = "a",
        Value = 1,
        Type = 'abs'
    )
    msg = "the following variables were used but not found in the model"
    msg2 = "regular expressions did not match any state variables or parameters to sum."
    expect_error(
        flexmodel(
            params = c(b = 1),
            state = c(X = 0),
            start_date = "2000-01-01",
            end_date = "2000-03-01",
            params_timevar = tv
        ),
        regexp = msg
    )
    expect_error(
        (
            flexmodel(
                params = c(a = 1),
                state = c(X = 0, Y = 0),
                start_date = "2000-01-01",
                end_date = "2000-01-01"
            )
            %>% add_rate("X", "Y", ~ (b) + (c))
        ),
        regexp = msg
    )
    expect_error(
        (
            flexmodel(
                params = c(a = 1),
                state = c(X = 0, Y = 0),
                start_date = "2000-01-01",
                end_date = "2000-01-01"
            )
            %>% add_state_param_sum("Zsum", "Z")
        ),
        regexp = msg2
    )
})

test_that("invalid state and parameter sum specification returns informative error msg", {
    reset_spec_version()
    tmb_mode()
    msg1 = "sums cannot have the same name as state variables or parameters"
    msg2 = "sum_name must be character-valued"
    msg3 = "can only specify one sum at a time"
    expect_error(
        (
            flexmodel(
                params = c(a = 1),
                state = c(X = 0, Y = 0),
                start_date = "2000-01-01",
                end_date = "2000-01-01"
            )
            %>% add_state_param_sum("X", "X")
        ),
        regexp = msg1
    )
    expect_error(
        (
            flexmodel(
                params = c(a = 1),
                state = c(X = 0, Y = 0),
                start_date = "2000-01-01",
                end_date = "2000-01-01"
            )
            %>% add_state_param_sum(1, "X")
        ),
        regexp = msg2
    )
    expect_error(
        (
            flexmodel(
                params = c(a = 1),
                state = c(X = 0, Y = 0),
                start_date = "2000-01-01",
                end_date = "2000-01-01"
            )
            %>% add_state_param_sum(c("X", "Y"), "X")
        ),
        regexp = msg3
    )
})

test_that("informative error is thrown if no rates are specified", {
    reset_spec_version()
    tmb_mode()
    msg = "no rates have been added to this model"
    expect_error(
        (
            flexmodel(
                params = c(a = 1),
                state = c(X = 0),
                start_date = "2000-01-01",
                end_date = "2000-01-02"
            )
            %>% update_tmb_indices()
        ),
        regexp = msg
    )
})

test_that("informative error is thrown if no rates are specified", {
    reset_spec_version()
    tmb_mode()
    op1 = options(MP_auto_outflow = FALSE)
    op2 = options(MP_warn_no_outflow = TRUE)
    msg = "model does not contain any outflow"
    expect_warning(
        (
            flexmodel(
                params = c(a = 1),
                state = c(X = 10, Y = 0),
                start_date = "2000-01-01",
                end_date = "2000-01-02"
            )
            %>% add_rate("X", "Y", ~ (a))
            %>% update_tmb_indices
        ),
        regexp = msg
    )
    options(op1)
    options(op2)
})

test_that("start_date <= end_date in flex models", {
  expect_error(
    flexmodel(
      params = read_params("ICU1.csv"),
      start_date = "2000-01-02",
      end_date = "2000-01-01"
    ),
    regexp = "start_date must be less than or equal to end_date"
  )
})

test_that("an error is thrown when params is not params_pansim and state is not provided", {
  expect_error(
    flexmodel(
      params = c(S = 0),
      start_date = "2000-01-01",
      end_date = "2000-01-02"
    ),
    regexp = "an initial state vector is required, because"
  )
})

test_that("pre-defined factors give the same answer as defining the rate with raw factors", {

  mm = flexmodel(
    params = c(a = 0.5, b = 0.25, c = 0.1),
    state = c(X = 1, Y = 2),
    start_date = "2000-01-01",
    end_date = "2000-02-05",
    do_make_state = FALSE
  )

  mm1 = (mm
    %>% add_rate("X", "Y", ~ (a) * (c) * (X) + (c) * (b) * (Y))
    %>% add_outflow %>% update_tmb_indices
  )

  mm2 = (mm
    %>% add_factr("alpha", ~ (a) * (X) + (b) * (Y))
    %>% add_rate("X", "Y", ~ (alpha) * (c))
    %>% add_outflow %>% update_tmb_indices
  )

  expect_equal(
    simulate_state_vector(mm1),
    simulate_state_vector(mm2)
  )

  mm = flexmodel(
      params = c(beta = 0.5, N = 100),
      state = c(S = 99, I = 1),
      start_date = "2000-01-01",
      end_date = "2000-01-05",
      do_make_state = FALSE
    )

  mm1 = (mm
    %>% add_factr("foi", ~ (beta) * (1/N) * (I))
    %>% add_rate("S", "I", ~ (foi))
    %>% add_outflow
    %>% update_tmb_indices
  )

  mm2 = (mm
    %>% add_rate("S", "I", ~ (beta) * (1/N) * (I))
    %>% add_outflow
    %>% update_tmb_indices
  )

  expect_equal(
    simulate_state_vector(mm1),
    simulate_state_vector(mm2)
  )

})

test_that("vector-valued pre-defined factors give consistent results", {
  strains = c("wild", "variant")
  state = c(
    S = 20000,
    I_wild = 49, I_variant = 1,
    R_wild = 0,   R_variant = 0
  )
  two_strain_model =
    flexmodel(
      params = c(
        gamma = 0.06,
        beta_wild = 0.15,
        beta_variant = 0.25,
        N = sum(state)
      ),
      state = state,
      start_date = "2000-01-01",
      end_date = "2000-01-02",
      do_hazard = TRUE,
      do_make_state = FALSE
    )

  two_strains_factr = (two_strain_model
    %>% vec_factr(
      "foi" %_% strains,
      vec("beta" %_% strains) * struc("1/N") * vec("I" %_% strains))
    %>% vec_rate("S", "I" %_% strains, vec("foi" %_% strains))
    %>% rep_rate("I", "R", ~ (gamma))
    %>% add_outflow()
    %>% update_tmb_indices
  )

  two_strains_no_factr = (two_strain_model
    %>% vec_rate(
      "S",
      "I" %_% strains,
      vec("beta" %_% strains) * struc("1/N") * vec("I" %_% strains)
    )
    %>% rep_rate("I", "R", ~ (gamma))
    %>% add_outflow()
    %>% update_tmb_indices
  )


  sims_factrs = simulate_state_vector(two_strains_factr)
  sims_no_factrs = simulate_state_vector(two_strains_no_factr)
  expect_equal(
    sims_factrs,
    sims_no_factrs
  )
})


test_that("sim_report expressions give correct results", {
  strains = c("wild", "variant")
  state = c(
    S = 20000,
    I_wild = 49, I_variant = 1,
    R_wild = 0,   R_variant = 0
  )
  two_strain_model =
    flexmodel(
      params = c(
        gamma = 0.06,
        beta_wild = 0.15,
        beta_variant = 0.25,
        N = sum(state),
        c_prop = 1e-1,
        c_delay_cv = 2.5e-1,
        c_delay_mean = 1.1
      ),
      state = state,
      start_date = "2000-01-01",
      end_date = "2000-01-02",
      do_hazard = TRUE,
      do_make_state = FALSE
    )

  model = (two_strain_model
     %>% vec_factr(
      "foi" %_% strains,
      vec("beta" %_% strains) * struc("1/N") * vec("I" %_% strains)
     )
     %>% vec_rate("S", "I" %_% strains, vec("foi" %_% strains))
     %>% rep_rate("I", "R", ~ (gamma))
     %>% add_sim_report_expr('report', ~ (I_wild) + (I_variant))
     %>% add_sim_report_expr('recov', ~ (S) * (S_to_I_wild) + (S) * (S_to_I_variant))
     %>% add_lag_diff("^report$")
     %>% add_conv("^recov$")
  )
  simulation_history(model)
})

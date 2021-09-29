library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)

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
spec_version <- "0.0.1"
print(spec_version)
options(MP_flex_spec_version = spec_version)
test_files <- "../../inst/tmb/"

cpp <- file.path(test_files, spec_version, "macpan.cpp")
dll <- file_path_sans_ext(cpp)
options(MP_flex_spec_dll = basename(dll))

compile(cpp)
dyn.load(dynlib(dll))

params <- read_params("ICU1.csv")
state <- McMasterPandemic::make_state(params = params)
M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

params <- read_params("ICU1.csv")
state <- make_state(params = params)
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

test_that("rate matrix types match", {
    expect_equal(
        class(tmb_sparse_ratemat),
        class(r_sparse_ratemat)
    )
})

test_that("rate matrix dimensions match", {
    expect_equal(
        tmb_sparse_ratemat@Dim,
        r_sparse_ratemat@Dim
    )
})

test_that("matrix elements are equal", {
    expect_equal(
        tmb_dense_ratemat,
        r_dense_ratemat
    )
})

dyn.unload(dynlib(dll))



spec_version <- "0.0.2"
print(spec_version)
options(MP_flex_spec_version = spec_version)
test_files <- "../../inst/tmb/"

cpp <- file.path(test_files, spec_version, "macpan.cpp")
dll <- file_path_sans_ext(cpp)
options(MP_flex_spec_dll = basename(dll))

compile(cpp)
dyn.load(dynlib(dll))


params <- read_params("ICU1.csv")
state <- make_state(params = params)
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
        %>% add_parallel_accumulators(c("X", "N", "P", "V"))
        %>% add_tmb_indices()
)

dd <- tmb_fun(test_model)

tmb_traj <- (state
    %>% c(dd$report()$concatenated_state_vector)
    %>% matrix(length(state), 4, dimnames = list(names(state), 1:4))
    %>% t()
    %>% as.data.frame()
)
r_traj <- run_sim_range(
    params = params, state = state, nt = 4,
    step_args = list(do_hazard = FALSE)
)[, names(state)]

test_that("simulated state trajectories are equal", {
    expect_equal(tmb_traj, r_traj)
})

dyn.unload(dynlib(dll))


spec_version <- "0.0.4"
print(spec_version)
options(MP_flex_spec_version = spec_version)

test_files <- "../../inst/tmb/"

cpp <- file.path(test_files, spec_version, "macpan.cpp")
dll <- file_path_sans_ext(cpp)
options(MP_flex_spec_dll = basename(dll))

compile(cpp)
dyn.load(dynlib(dll))

params <- read_params("ICU1.csv")
state <- make_state(params = params)
M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

tv_dat <- data.frame(
    Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.5, 0.1, 0.05),
    Type = c("rel_prev", "rel_orig", "rel_prev")
)

test_model <- (init_model(
    params, state,
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

tmb_traj <- (state
    %>% c(dd$report()$concatenated_state_vector)
    %>% matrix(length(state), test_model$iters + 1,
        dimnames = list(names(state), 1:(test_model$iters + 1))
    )
    %>% t()
    %>% as.data.frame()
)

r_traj <- run_sim(
    params = params, state = state,
    start_date = test_model$start_date,
    end_date = test_model$end_date,
    params_timevar = tv_dat,
    condense = FALSE,
    step_args = list(do_hazard = FALSE)
)[, names(state)] %>%
    as.data.frame()

row.names(r_traj) <- row.names(tmb_traj) <- NULL

test_that("simulated state trajectories are equal", {
    expect_equal(tmb_traj, r_traj)
})




spec_version <- "0.0.5"
print(spec_version)
options(MP_flex_spec_version = spec_version)

test_files <- "../../inst/tmb/"

cpp <- file.path(test_files, spec_version, "macpan.cpp")
dll <- file_path_sans_ext(cpp)
options(MP_flex_spec_dll = basename(dll))

compile(cpp)
dyn.load(dynlib(dll))

params <- read_params("ICU1.csv")
state <- make_state(params = params)
M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

tv_dat <- data.frame(
    Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.5, 0.1, 0.05),
    Type = c("rel_prev", "rel_orig", "rel_prev")
)


test_model <- (init_model(
    params, state,
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
    %>% matrix(length(state), test_model$iters + 1,
        dimnames = list(names(state), 1:(test_model$iters + 1))
    )
    %>% t()
    %>% as.data.frame()
)

r_traj <- run_sim(
    params = params, state = state,
    start_date = test_model$start_date,
    end_date = test_model$end_date,
    params_timevar = tv_dat,
    condense = FALSE,
    step_args = list(do_hazard = TRUE)
)[, names(state)] %>%
    as.data.frame()

row.names(r_traj) <- row.names(tmb_traj) <- NULL

test_that("simulated state trajectories are equal", {
    expect_equal(tmb_traj, r_traj)
})

test_that("mock objective function is correct", {
    expect_equal(sum(tmb_traj[31, ]), dd$fn(dd$par))
})

test_that("tmb-computed gradient equals numerical gradient", {

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
})


test_that("use_flex flag does not change results", {
    tmb_sim <- run_sim(
        params = params, state = state,
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

    # exceptions to drop-in replacement:
    # 1. don't require that the attributes are in the same order
    # 2. don't require that the r version returns everything that the tmb
    #    version does (e.g. flexmodel)
    # 3. don't require that the row.names are identical (is this ok?
    #    the r version counts iterations with skips, but is this informative?)
    # 4. don't require that the call is identical (obvious i guess, but
    #    being exhaustive)
    attr(tmb_sim, 'row.names') = attr(r_sim, 'row.names') = NULL
    attr(tmb_sim, 'call') = attr(r_sim, 'call') = NULL
    for(a in names(attributes(r_sim))) {
        expect_equal(attr(tmb_sim, a), attr(r_sim, a))
    }

    attributes(r_sim) <- attributes(tmb_sim) <- NULL
    expect_equal(tmb_sim, r_sim)
})



spec_version <- "0.0.6"
print(spec_version)
options(MP_flex_spec_version = spec_version)

test_files <- "../../inst/tmb/"

cpp <- file.path(test_files, spec_version, "macpan.cpp")
dll <- file_path_sans_ext(cpp)
options(MP_flex_spec_dll = basename(dll))

compile(cpp)
dyn.load(dynlib(dll))

test_that('time-varying parameters are correctly updated on C++ side', {
    params <- read_params("ICU1.csv")
    state <- make_state(params = params)
    M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

    tv_dat <- data.frame(
        Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
        Symbol = c("beta0", "beta0", "beta0"),
        Value = c(0.5, 0.1, 0.05),
        Type = c("rel_prev", "rel_orig", "rel_prev")
    )

    mm = make_unflexmodel(
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
        params = params, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        params_timevar = tv_dat,
        condense = FALSE,
        step_args = list(do_hazard = TRUE),
        use_flex = FALSE
    )

    attr(tmb_sim, 'row.names') = attr(r_sim, 'row.names') = NULL
    attr(tmb_sim, 'call') = attr(r_sim, 'call') = NULL
    for(a in names(attributes(r_sim))) {
        print(a)
        expect_equal(attr(tmb_sim, a), attr(r_sim, a))
    }

    attributes(r_sim) <- attributes(tmb_sim) <- NULL
    expect_equal(tmb_sim, r_sim)
})

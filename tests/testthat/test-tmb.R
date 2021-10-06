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

    compare_sims(r_sim, tmb_sim)
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

test_that('it remains ok to _not_ use time-varying parameters', {
    params <- read_params("ICU1.csv")
    state <- make_state(params = params)
    M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

    tmb_sim <- run_sim(
        params = params, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        condense = FALSE,
        step_args = list(do_hazard = TRUE),
        use_flex = TRUE
    )

    r_sim <- run_sim(
        params = params, state = state,
        start_date = "2021-09-10",
        end_date = "2021-10-10",
        condense = FALSE,
        step_args = list(do_hazard = TRUE),
        use_flex = FALSE
    )
    compare_sims(r_sim, tmb_sim)
})




spec_version <- "0.1.0"
print(spec_version)
options(MP_flex_spec_version = spec_version)

test_files <- "../../inst/tmb/"

# TODO: not yet using spec-version-specific dll because it is not yet created
#cpp <- file.path(test_files, spec_version, "macpan.cpp")
#dll <- file_path_sans_ext(cpp)
#options(MP_flex_spec_dll = basename(dll))
options(MP_flex_spec_dll = "McMasterPandemic")

#compile(cpp)
#dyn.load(dynlib(dll))

test_that('simple models still work when structure is allowed', {
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

test_that("simple vaccination model in TMB matches and is faster than existing R model", {
    options(macpan_pfun_method = "grep")
    params <- read_params("ICU1.csv")
    state <- make_state(params = params)
    vax_params <- expand_params_vax(
        params = params,
        model_type = "twodose"
    )
    vax_state <- expand_state_vax(
        x = state,
        model_type = "twodose",
        unif = FALSE
    )

    # problem dimensions
    (epi_states = c(attr(vax_state, "epi_cat")))
    (asymp_cat = c("S", "E", "Ia", "Ip", "R"))
    (vax_cat = c(attr(vax_state, "vax_cat")))
    (dose_from = rep(asymp_cat, 2))
    (dose_to = c(asymp_cat, rep("V", 5)))

    # Specify structure of the force of infection calculation
    Istate = (McMasterPandemic:::expand_names(
        c('Ia', 'Ip', 'Im', 'Is'),   # = Icats
        attr(vax_params, 'vax_cat')) # = vax_cats
        %>% struc
    )
    baseline_trans_rates =
        struc(
            'Ca',
            'Cp',
            '(1 - iso_m) * (Cm)',
            '(1 - iso_s) * (Cs)') *
        struc('(beta0) * (1/N)')
    vax_trans_red = struc_block(struc(
        '1',
        '1',
        '(1 - vax_efficacy_dose1)',
        '(1 - vax_efficacy_dose1)',
        '(1 - vax_efficacy_dose2)'),
        row_times = 1, col_times = 5)

    alpha   = c("alpha", "alpha", "vax_alpha_dose1", "vax_alpha_dose1", "vax_alpha_dose2")
    mu      = c("mu",    "mu",    "vax_mu_dose1",    "vax_mu_dose1",    "vax_mu_dose2")
    sigma   = struc("sigma")
    gamma_p = struc("gamma_p")
    E_to_Ia_rates  = struc(           alpha ) * sigma
    E_to_Ip_rates  = struc(complement(alpha)) * sigma
    Ip_to_Im_rates = struc(              mu ) * gamma_p
    Ip_to_Is_rates = struc(complement(   mu)) * gamma_p

    test_model <- (init_model(
        vax_params, vax_state,
        start_date = "2021-09-09", end_date = "2021-10-09"
    )

    # Flow within vaccination categories,
    # with constant rates across categories
    %>% rep_rate("Ia",   "R",    ~                      (gamma_a))
    %>% rep_rate("Im",   "R",    ~                      (gamma_m))
    %>% rep_rate("Is",   "D",    ~ (    nonhosp_mort) * (gamma_s))
    %>% rep_rate("Is",   "H",    ~ (1 - nonhosp_mort) * (gamma_s) * (    phi1))
    %>% rep_rate("Is",   "X",    ~ (1 - nonhosp_mort) * (gamma_s) * (    phi1))
    %>% rep_rate("Is",   "ICUs", ~ (1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (1 - phi2))
    %>% rep_rate("Is",   "ICUd", ~ (1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (    phi2))
    %>% rep_rate("ICUs", "H2",   ~                                  (    psi1))
    %>% rep_rate("ICUd", "D",    ~                                  (    psi2))
    %>% rep_rate("H2",   "R",    ~                                  (    psi3))
    %>% rep_rate("H",    "R",    ~ (rho))

    # Flow within vaccination categories,
    # with rates that depend on category
    # (see struc objects created above)
    %>% vec_rate("E", "Ia",  E_to_Ia_rates)
    %>% vec_rate("E", "Ip",  E_to_Ip_rates)
    %>% vec_rate("Ip", "Im", Ip_to_Im_rates)
    %>% vec_rate("Ip", "Is", Ip_to_Is_rates)

    # Vaccination Response Rates
    %>% add_rate("R_vaxdose1", "R_vaxprotect1",  ~ (vax_response_rate_R))
    %>% add_rate("R_vaxdose2", "R_vaxprotect2",  ~ (vax_response_rate_R))
    %>% add_rate("S_vaxdose1", "S_vaxprotect1",  ~ (vax_response_rate))
    %>% add_rate("S_vaxdose2", "S_vaxprotect2",  ~ (vax_response_rate))

    # Forces of Infection
    %>% vec_rate(
        "S" %_% vax_cat,
        "E" %_% vax_cat,
        kronecker(vax_trans_red, t(baseline_trans_rates)) %*% Istate
    )

    # Sums across vaccination categories
    %>% add_state_param_sum("asymp_unvax_N",       asymp_cat %_% "unvax")
    %>% add_state_param_sum("asymp_vaxprotect1_N", asymp_cat %_% "vaxprotect1")

    # Flow among vaccination categories
    # (see dose_* above for epi states that are involved)
    %>% rep_rate(
        dose_from %_% 'unvax',
        dose_to   %_% 'vaxdose1',
        ~ (    vax_prop_first_dose) * (vax_doses_per_day) * (1 / asymp_unvax_N))
    %>% rep_rate(
        dose_from %_% 'vaxprotect1',
        dose_to   %_% 'vaxdose2',
        ~ (1 - vax_prop_first_dose) * (vax_doses_per_day) * (1 / asymp_vaxprotect1_N))

    # Technical steps
    %>% add_parallel_accumulators(c('V' %_% vax_cat, 'X' %_% vax_cat))
    %>% add_tmb_indices()
    )
    tmb_strt = Sys.time()
    tmb_sim <- run_sim(
        params = vax_params, state = vax_state,
        start_date = "2021-09-09", end_date = "2021-10-09",
        condense = FALSE,
        step_args = list(do_hazard = TRUE),
        flexmodel = test_model,
        use_flex = TRUE
    )
    tmb_nd = Sys.time()
    r_strt = Sys.time()
    r_sim = run_sim(
        params = vax_params, state = vax_state,
        start_date = "2021-09-09", end_date = "2021-10-09",
        condense = FALSE,
        step_args = list(do_hazard = TRUE)
    )
    r_nd = Sys.time()
    tmb_speed = as.numeric(tmb_nd - tmb_strt)
    r_speed = as.numeric(r_nd - r_strt)
    print(paste0('tmb speed-up: ', round(r_speed/tmb_speed), 'x'))
    expect_lt(tmb_speed, r_speed)
    compare_sims(r_sim, tmb_sim)
})

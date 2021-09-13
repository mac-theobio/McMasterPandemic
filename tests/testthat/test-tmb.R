library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)

# TODO/FIXME/BACKGROUND --------------------------------------------
# rename/repurpose this once the tmb engine is properly added
# to the package. for now refer to a particular version
# of the spec and store the C++ files in the
# appropriate directory within inst. also for now, assume that this
# test file is being run inside the directory in which it should
# live within the package structure
# (i.e. inside McMasterPandemic/tests/testthat). this should
# also work for `make pkgtest` because tests are run in the
# directory of the test file (as far as I can tell), but it could
# be fragile in other contexts? in any case, once we have
# the tmb-refactored code in the package, we will use the
# installed dlls and R-side utilities but retain the structure of
# these tests by spec version. possibly we can roll-up all up all
# spec versions with the same major version together (i.e. those
# that are backwards compatible).
# see https://canmod.net/misc/flex_specs for more on
# spec versioning.
# ------------------------------------------------------------------
spec_version = "0.0.1"
options(MP_flex_spec_version = spec_version)
test_files = "../../inst/tmb/"

cpp = file.path(test_files, spec_version, "macpan.cpp")
dll = file_path_sans_ext(cpp)

compile(cpp)
dyn.load(dynlib(dll))

params = read_params('ICU1.csv')
state = McMasterPandemic::make_state(params = params)
M = McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

params = read_params('ICU1.csv')
state = make_state(params = params)
test_model = (
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
  %>% add_rate("S",  "E", ~
                 (Ia) * (beta0) * (1/N) * (Ca) +
                 (Ip) * (beta0) * (1/N) * (Cp) +
                 (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
                 (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
  %>% add_tmb_indices()
)

indices = test_model$tmb_indices$make_ratemat_indices
dd <- MakeADFun(data = list(state = c(test_model$state),
                            ratemat = M,
                            from = indices$from,
                            to = indices$to,
                            count = indices$count,
                            spi = indices$spi,
                            modifier = indices$modifier),
                parameters = list(params=c(test_model$params)),
                DLL=basename(dll))

tmb_sparse_ratemat = dd$report()$ratemat
r_dense_ratemat = as.matrix(M)
tmb_dense_ratemat = as.matrix(tmb_sparse_ratemat)
r_sparse_ratemat = M
dimnames(tmb_sparse_ratemat) <-
  dimnames(r_sparse_ratemat) <-
  dimnames(tmb_dense_ratemat) <-
  dimnames(r_dense_ratemat) <- NULL

test_that("rate matrix types match", {
  expect_equal(
    class(tmb_sparse_ratemat),
    class(r_sparse_ratemat))
})

test_that("rate matrix dimensions match", {
  expect_equal(
    tmb_sparse_ratemat@Dim,
    r_sparse_ratemat@Dim)
})

test_that("matrix elements are equal", {
  expect_equal(
    tmb_dense_ratemat,
    r_dense_ratemat)
})





spec_version = "0.0.2"
options(MP_flex_spec_version = spec_version)
test_files = "../../inst/tmb/"

cpp = file.path(test_files, spec_version, "macpan.cpp")
dll = file_path_sans_ext(cpp)

compile(cpp)
dyn.load(dynlib(dll))


params = read_params('ICU1.csv')
state = make_state(params = params)
M = McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

test_model = (
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
  %>% add_rate("S",  "E", ~
                 (Ia) * (beta0) * (1/N) * (Ca) +
                 (Ip) * (beta0) * (1/N) * (Cp) +
                 (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
                 (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
  %>% add_parallel_accumulators(c("X", "N", "P", "V"))
  %>% add_tmb_indices()
)

indices = test_model$tmb_indices$make_ratemat_indices
update_indices = test_model$tmb_indices$update_ratemat_indices
par_accum_indices = test_model$tmb_indices$par_accum_indices

numIters = 3 # This value should have been kept in test_model

dd <- MakeADFun(data = list(state = c(test_model$state),
                            ratemat = M,
                            from = indices$from,
                            to = indices$to,
                            count = indices$count,
                            spi = indices$spi,
                            modifier = indices$modifier,
                            update_from = update_indices$from,
                            update_to = update_indices$to,
                            update_count = update_indices$count,
                            update_spi = update_indices$spi,
                            update_modifier = update_indices$modifier,
                            par_accum_indices = par_accum_indices,
                            numIterations = numIters),
                parameters = list(params=c(test_model$params)),
                DLL=basename(dll))

tmb_traj = (state
  %>% c(dd$report()$concatenated_state_vector)
  %>% matrix(length(state), 4, dimnames = list(names(state), 1:4))
  %>% t
  %>% as.data.frame
)
r_traj = run_sim_range(
  params = params, state = state, nt = 4,
  step_args = list(do_hazard = FALSE))[,names(state)]

test_that("simulated state trajectories are equal", {
  expect_equal(tmb_traj, r_traj)
})




spec_version = "0.0.4"
options(MP_flex_spec_version = spec_version)

test_files = "../../inst/tmb/"

cpp = file.path(test_files, spec_version, "macpan.cpp")
dll = file_path_sans_ext(cpp)

compile(cpp)
dyn.load(dynlib(dll))

params = read_params('ICU1.csv')
state = make_state(params = params)
M = McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

tv_dat = data.frame(
  Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
  Symbol = c("beta0", "beta0", "beta0"),
  Value = c(0.5, 0.1, 0.05),
  Type = c('rel_prev', 'rel_orig', 'rel_prev')
)

test_model = (init_model(
  params, state,
  start_date = "2021-09-09", end_date = "2021-10-09",
  timevar_piece_wise = tv_dat)
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
  %>% add_rate("S",  "E", ~
                 (Ia) * (beta0) * (1/N) * (Ca) +
                 (Ip) * (beta0) * (1/N) * (Cp) +
                 (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
                 (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
  %>% add_parallel_accumulators(c("X", "N", "P", "V"))
  %>% add_tmb_indices()
)

from = test_model$tmb_indices$make_ratemat_indices$from
to = test_model$tmb_indices$make_ratemat_indices$to
count = test_model$tmb_indices$make_ratemat_indices$count
spi = test_model$tmb_indices$make_ratemat_indices$spi
modifier = test_model$tmb_indices$make_ratemat_indices$modifier

updateidx = test_model$tmb_indices$updateidx
breaks = test_model$timevar$piece_wise$breaks
count_of_tv_at_breaks = test_model$timevar$piece_wise$count_of_tv_at_breaks
tv_spi = test_model$timevar$piece_wise$schedule$tv_spi
tv_val = test_model$timevar$piece_wise$schedule$tv_val
par_accum_indices = test_model$tmb_indices$par_accum_indices
iters = test_model$iters

dd <- MakeADFun(data = list(state = c(state),
                            ratemat = M,
                            from = from,
                            to = to,
                            count = count,
                            spi = spi,
                            modifier = modifier,
                            updateidx = c(updateidx),
                            breaks = breaks,
                            count_of_tv_at_breaks = count_of_tv_at_breaks,
                            tv_spi = tv_spi,
                            tv_val = tv_val,
                            par_accum_indices = par_accum_indices,
                            numIterations = iters),
                parameters = list(params=c(test_model$params)),
                DLL=basename(dll))

tmb_traj = (state
            %>% c(dd$report()$concatenated_state_vector)
            %>% matrix(length(state), test_model$iters + 1,
                       dimnames = list(names(state), 1:(test_model$iters+1)))
            %>% t
            %>% as.data.frame
)

r_traj = run_sim(
  params = params, state = state,
  start_date = test_model$start_date,
  end_date = test_model$end_date,
  params_timevar = tv_dat,
  condense = FALSE,
  step_args = list(do_hazard = FALSE))[,names(state)] %>%
  as.data.frame

row.names(r_traj) = row.names(tmb_traj) = NULL

test_that("simulated state trajectories are equal", {
  expect_equal(tmb_traj, r_traj)
})




spec_version = "0.0.5"
options(MP_flex_spec_version = spec_version)

test_files = "../../inst/tmb/"

cpp = file.path(test_files, spec_version, "macpan.cpp")
dll = file_path_sans_ext(cpp)

compile(cpp)
dyn.load(dynlib(dll))

params = read_params('ICU1.csv')
state = make_state(params = params)
M = McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

tv_dat = data.frame(
  Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
  Symbol = c("beta0", "beta0", "beta0"),
  Value = c(0.5, 0.1, 0.05),
  Type = c('rel_prev', 'rel_orig', 'rel_prev')
)

test_model = (init_model(
  params, state,
  start_date = "2021-09-09", end_date = "2021-10-09",
  timevar_piece_wise = tv_dat)
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
  %>% add_rate("S",  "E", ~
                 (Ia) * (beta0) * (1/N) * (Ca) +
                 (Ip) * (beta0) * (1/N) * (Cp) +
                 (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
                 (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
  %>% add_parallel_accumulators(c("X", "N", "P", "V"))
  %>% add_tmb_indices()
)

#sp = c(test_model$state, test_model$params)
from = test_model$tmb_indices$make_ratemat_indices$from
to = test_model$tmb_indices$make_ratemat_indices$to
count = test_model$tmb_indices$make_ratemat_indices$count
spi = test_model$tmb_indices$make_ratemat_indices$spi
modifier = test_model$tmb_indices$make_ratemat_indices$modifier

updateidx = test_model$tmb_indices$updateidx
breaks = test_model$timevar$piece_wise$breaks
count_of_tv_at_breaks = test_model$timevar$piece_wise$count_of_tv_at_breaks
tv_spi = test_model$timevar$piece_wise$schedule$tv_spi
tv_val = test_model$timevar$piece_wise$schedule$tv_val
par_accum_indices = test_model$tmb_indices$par_accum_indices

do_hazard = TRUE

dd <- MakeADFun(data = list(state = c(state), #c(test_model$state),
                            ratemat = M,
                            from = from,
                            to = to,
                            count = count,
                            spi = spi,
                            modifier = modifier,
                            updateidx = c(updateidx),
                            breaks = breaks,
                            count_of_tv_at_breaks = count_of_tv_at_breaks,
                            tv_spi = tv_spi,
                            tv_val = tv_val,
                            par_accum_indices = par_accum_indices,
                            do_hazard = do_hazard,
                            numIterations = test_model$iters),
                parameters = list(params=c(test_model$params)),
                DLL=basename(dll))


tmb_traj = (state
            %>% c(dd$report()$concatenated_state_vector)
            %>% matrix(length(state), test_model$iters + 1,
                       dimnames = list(names(state), 1:(test_model$iters+1)))
            %>% t
            %>% as.data.frame
)

r_traj = run_sim(
  params = params, state = state,
  start_date = test_model$start_date,
  end_date = test_model$end_date,
  params_timevar = tv_dat,
  condense = FALSE,
  step_args = list(do_hazard = do_hazard))[,names(state)] %>%
  as.data.frame

row.names(r_traj) = row.names(tmb_traj) = NULL

test_that("simulated state trajectories are equal", {
  expect_equal(tmb_traj, r_traj)
})

test_that("mock objective function is correct", {
  expect_equal(sum(tmb_traj[31,]), dd$fn(dd$par))
})

test_that("tmb-computed gradient equals (with tolerance) a numerically-approximated gradient", {
  numeric_gradient = numDeriv::grad(dd$fn, dd$par)
  tmb_gradient = dd$gr(dd$par)
  attributes(numeric_gradient) = attributes(tmb_gradient) = NULL
  expect_equal(tmb_gradient, numeric_gradient,
               tolerance = 1e-04) # test fails with tol=1e-05
})

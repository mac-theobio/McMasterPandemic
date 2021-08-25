library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)

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
context(paste0("tests for spec version: ", spec_version))
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
context(paste0("tests for spec version: ", spec_version))
test_files = "../../inst/tmb/"

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
  %>% add_parallel_accumulators(c("X", "N", "P", "V"))
  %>% add_tmb_indices()
)

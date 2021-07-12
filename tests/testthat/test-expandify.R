test_that("states are properly expanded with vaccination", {
    params <- read_params("PHAC_testify.csv")
    ss <- make_state(params = params, testify = FALSE)
    ssv <- expand_state_vax(ss)
    ## check expand
    expect_equal(sum(ssv), sum(ss))
    ## check condense
    ## expect_equal(condense_vax(ssv), ss, ignore_attr = TRUE)
})

## test_that("states are properly expanded with age *and* vaccination", {
##   age_cat <- mk_agecats(0, 80, da = 30)
##   Nvec <- rep(1e5, length(age_cat))
##   base_params <- update(
##     read_params("PHAC_testify.csv")
##     , N = sum(Nvec)
##     , testing_intensity = 0
##   )
##   age_params <- expand_params_age(
##     params = base_params,
##     age_cat = age_cat,
##     Nvec = Nvec
##   )
##   age_vax_params <- expand_params_vax(age_params)
##
##   make_state(params = age_vax_params)
## })

## TEST IDEAS

## make ratemat for just vax case, and check that on the 3 diagonal blocks (corresponding to the three vax statuses), where the top two blocks are exactly equal and the bottom one differs only in vaccine efficacy and whatever changes we make to the epidemiological parameters for vaccinated individuals

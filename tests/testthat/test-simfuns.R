library(testthat)
library(McMasterPandemic)

n <- 15
names <- c(99:(99 + n - 1))

M <- matrix(runif(n * n), nrow = n, ncol = n)
rownames(M) <- names
colnames(M) <- names
v <- c(1:n)
names(v) <- names


test_that("New sweep works correctly", {
    baseline <- sweep(M, v, MARGIN = 1, FUN = "*")
    expect(identical(col_multiply(M, v), baseline), "col_multiply not working as expected!")
})

params1 <- read_params("ICU1.csv")
state1 <- make_state(params = params1)
sdate <- "2020-02-10"
edate <- "2020-06-01"

test_that("State variable <0 warning doesn't throw in the typical case", {
  expect_warning(run_sim(params = params1, state = state1, start_date = sdate, end_date = edate,
                         step_args = list(do_hazard = FALSE)), regexp = NA)
})

state1["S"] <- 0
params1["N"] <- sum(state1)

test_that("State variable <0 warning works correctly", {
    ## Note: if the code is changed to prevent things from going below 0, this test will fail.
    expect_warning(run_sim(params = params1, state = state1, start_date = sdate, end_date = edate,
                          step_args = list(do_hazard = FALSE))),
            label = "Note: if the code is changed to prevent things from going below 0, this test will fail. Change the test if this happens.", regexp = "One or more state variables is negative")
})


test_that(" No warnings thrown if state variables are all capital letters", {
    expect_warning(state1 <- make_state(params = params1), regexp = NA, label = "Inappropriate warning for non-capitalized state when state variables are capitalized correctly")
})

test_that("Warnings thrown if state variables are not all capital letters", {
    expect_warning(state1 <- make_state(params = params1, type = "test_warning_throw",
                                       use_eigvec = FALSE),
                   regexp = "Not all state variables are capital letters", label = "Non-capitalized state warning not being thrown")
})


## microbenchmark(baseline = sweep(M, v, MARGIN=1, FUN="*"))
## microbenchmark({new =t(t(M) %*% diag(v));
## n = names(v);
## rownames(new) = n})
## microbenchmark(new= t(t(M) * rep(v, rep.int(nrow(M), length(v)))))

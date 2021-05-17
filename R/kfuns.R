## kernel functions
## not exported or documented yet!

## FIXME: how to aggregate/thin/adjust time scale?
## compute kernel by brute force
transKernel <- function(par, steps = 100, do_hazard = TRUE,
                        ndt = 1) {
    if (ndt > 1) warning("ndt not fully implemented")
    par[["N"]] <- 1 ## ? redundant ?
    state <- make_state(N = 1, E0 = 1, use_eigvec = FALSE)
    return(run_sim_range(par, state,
        nt = steps * ndt,
        step_args = list(do_hazard = do_hazard)
    ))
}
## allow caching of results
if (requireNamespace("memoise")) {
    transKernel <- memoise::memoise(transKernel)
}

## FIXME: kernel should ideally be an object with k and lag
## (a class with methods)
## Building with r instead of <U+03BB> for this reason

## Investigate r (combine with uniroot to get Euler's r)
discountGap <- function(r, k) {
    lag <- 1:length(k)
    discountR <- sum(k * exp(-lag * r))
    return(discountR - 1)
}

## Investigate <U+03BA> (combine with uniroot to get <U+03BA>_eff)
kappaGap <- function(kappa, rho, R) {
    R_est <- (1 + rho * kappa)^(1 / kappa)
    return(R - R_est)
}

## moments of a kernel with parameters/convolution vector k
## FIXME: principled way to deal with lower-bound issues
kernelMoments <- function(k, lwr = 0.01) {
    lag <- seq(length(k))
    R0 <- sum(k)
    Gbar <- sum(k * lag) / R0
    Gvar <- sum(k * lag^2) / R0 - Gbar^2
    r0 <- (uniroot(discountGap, k = k, lower = lwr, upper = 2))$root
    kappa_eff <- (uniroot(kappaGap,
        rho = r0 * Gbar,
        R = R0, lower = lwr, upper = 2
    ))$root

    return(c(
        R0 = R0, Gbar = Gbar, r0 = r0,
        kappa = Gvar / Gbar^2,
        kappa_eff = kappa_eff
    ))
}

## What to multiply beta by to get a desired r
##' @importFrom stats uniroot
rmult <- function(k, r) {
    uniroot(
        f = function(m) {
            discountGap(r, m * k)
        },
        lower = 1 / 10, upper = 10
    )$root
}

##
##' run a pure-exponential sim;
##' uses run_sim_range with a population of 1 (proportions) and a very small starting value,
##'  run for 100 steps (by default)
##' used to calculate either r (technically r0) or eigenvector (for distributing initial exposed across states)
##' @param params parameters
##' @param steps number of steps to run
##' @param ndt sub-time steps
##' @param do_hazard run with hazard
##' @param testify testing compartments
##' @param return_val return growth rate or eigenvector?
##' @examples
##' pp <- read_params("PHAC_testify.csv")
##' rExp(pp)
##' rExp(pp,return_val="eigenvector")
##' rExp(pp,return_val="eigenvector",testify=TRUE)
##' @export
rExp <- function(params, steps = 100, ndt = 1,
                 do_hazard = FALSE,
                 testify = has_testing(params = params),
                 return_val = c("r0", "eigenvector", "sim")) {
    return_val <- match.arg(return_val)
    if (ndt > 1) warning("ndt not fully implemented")
    params[["N"]] <- 1 ## ? redundant ?
    ## potential recursion here: have to make sure
    state <- make_state(
        N = 1, E0 = 1e-5, type = "ICU1",
        use_eigvec = FALSE,
        params = params,
        testify = FALSE
    ) ## FIXME: don't assume ICU1?
    M <- make_ratemat(state = state, params = params)
    if (testify) {
        M <- testify(M, params)
        state <- expand_stateval_testing(state, method = "untested")
    }
    r <- run_sim_range(params,
        state,
        nt = steps * ndt,
        M = M,
        step_args = list(
            do_hazard = do_hazard,
            do_exponential = TRUE,
            testwt_scale = "none"
        )
    )
    if (return_val == "sim") {
        return(r)
    }
    nn <- ndt * steps
    ## DRY: get_evec()
    drop_vars <- c("date", "t", "D", "foi", "R", "X", "N", "P")
    if (!testify) drop_vars <- c("S", drop_vars)
    drop_re <- paste0("^(", paste(drop_vars, collapse = "|"), ")")
    ## FIXME: safer version of this:
    ##   keep_vars_regexp <- "^[EIHh]"
    ##   unlist(x[pos, grepl(keep_vars_regexp, names(r))])
    uf <- function(x, pos) unlist(x[pos, !grepl(drop_re, names(r))])
    r_last <- uf(r, nn)
    r_nextlast <- uf(r, nn - 1)
    ret <- switch(return_val,
        ## log mean(x(t+1)/x(t))
        r0 = mean(log(r_last / r_nextlast)),
        ## normalized state vector at last time step
        eigenvector = unlist(r_last / sum(r_last))
    )
    return(ret)
}

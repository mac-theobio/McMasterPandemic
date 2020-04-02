##' Back-calculate initial conditions
##' Suppose we have a slope and an intercept for some variable
##' (e.g. hospitalizations) starting
##' from day t1.  We want an initial condition for day t0.
##' (Don't go too far back or rounding will kill us)#
#' @param date0 date to get initial conditions for
##' @param date1 initial date for regressions (t=0)
##' @param params a parameter set
##' @param int regression intercept
##' @param slope regression slope
##' @param var regression variable
##' @examples get_init(var="H",int=2.85,slope=0.176, param=params)
##' @export
get_init <- function(date0=ldmy("1-Mar-2020"),
                     date1=ldmy("25-Mar-2020"),
                     params,
                     int,
                     slope,
                     var="H") {
    ldmy <- lubridate::dmy  ## sugar
    state <- make_state(N=params[["N"]],E0=0.001)
    J <- make_jac(state=state,params=params)    
    ee <- eigen(J)
    v <- ee$vectors
    rownames(v) <- names(state)
    (dom_vec <- v[,which.max(ee$values)])
    dom_vec <- abs(dom_vec[-1]) ## drop S
    ## back-project modeled var from date1 to date0
    ## (predicted log(modeled var) on date0
    int0 <- int - as.numeric(date1-date0)*slope
    expected <- round(exp(int0)*dom_vec/dom_vec[[var]])
    expected <- expected[expected>0]
    state[names(expected)] <- expected
    ## not important, but adjust to keep same total pop
    state[["S"]] <- params[["N"]] - sum(expected)
    return(state)
}



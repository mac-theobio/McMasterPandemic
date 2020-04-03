## FIXME: generalize to take vector of params (log change)
##  plus parallel list of parameter vectors e.g.
##  log_delta=c(beta0=...,lambda=...)
##  pars_adj= list("beta0",c("gamma","lambda_s","lambda_m", "lambda_a"))
adjust_params <- function(log_d_beta0,log_d_lambda,params) {
    lambda_vars <- c("gamma","lambda_s","lambda_m", "lambda_a")
    params[lambda_vars] <- params[lambda_vars]*exp(log_d_lambda)
    params[["beta0"]] <- params[["beta0"]]*exp(log_d_beta0)
    return(params)
}

## FIXME: generalize s_new based on names of target
##  (r, R0, Gbar, kappa)
badness <- function(delta, params, target) {
    p_new <- adjust_params(delta[1],delta[2],params)
    s_new <- rep(NA,length(target))
    names(s_new) <- names(target)
    GIm <- get_GI_moments(p_new)
    for (i in names(target)) {
        s_new[[i]] <- switch(i,
                              r=get_r(p_new),
                              R0=get_R0(p_new),
                              Gbar=GIm[["Gbar"]],
                              kappa=GIm[["kappa"]])
    }
    res <- sum((s_new-target)^2)
    ## cat(delta,s_new,target,res,"\n")
    return(res)
}

## round-trip: should be zero
## badness(c(0,0), params=p1,target=c(r=get_r(p1,state), R0=get_R0(p1)),
##        state=state)

## 
## badness(c(0,0), params=p1,target=c(r=0.23,R0=3), state=state)

##' adjust pars to match targets
##' @export
fix_pars <- function(params, target=c(r=0.23,R0=3)) {
    ## cc <- emdbook::curve3d(badness(c(x,y),params,target=target),
    ##                        xlim=c(-1,1),ylim=c(-1,1),
    ##                        sys3d="image")
    opt1 <- optim(par=c(0,0), fn= badness, method="Nelder-Mead",
                  target=target, params=params)
    ##                  control=list(trace=TRUE))
    p_new <- adjust_params(opt1$par[1],opt1$par[2], params)
    return(p_new)
}

##' given a log-linear intercept and slope and target values for R0 or GI,
##' return initial conditions and parameters that will (approximately?)
##' recapitulate the observed log-linear dynamics
##' @param int log-linear intercept (expected value on date1)
##' @param slope log-linear slope (growth rate)
##' @param target vector of target statistics (allowed: R0, Gbar (mean generation interval), or kappa (CV^2 of generation interval)
##' @param date0 date for initial conditions
##' @param date1 initial date for regression
##' @param var variable for regression (H, D, ICU)
##' @examples
##' params <- read_params(system.file("params","ICU1.csv",package="McMasterPandemic"))
##' cc1 <- calibrate(int=log(200),slope=0.23,pop=1e6,params)
##' summary(cc1$params)
##' cc2 <- calibrate(int=log(200),slope=0.23,pop=1e6,params,target=c(R0=3))
##' summary(cc2$params)
## FIXME: warn if starting conditions are too low 
##' 
##' @export
calibrate <- function(int,slope,pop,params,
                      target=c(Gbar=6),
                      date0="1-Mar-2020",
                      date1="25-Mar-2020",
                      var="H") {
    date0 <- ldmy(date0); date1 <- ldmy(date1)
    ok_targets <- c("Gbar","kappa","R0")
    bad_targets <- setdiff(names(target),ok_targets)
    if (length(bad_targets)>0) stop("bad targets: ",paste(bad_targets,collapse=", "))
    ## first adjust base params to set r and R0 (or GI) to specified value
    params[["N"]] <- pop
    state <- make_state(N=pop,E0=1)
    target[["r"]] <- slope
    p2 <- fix_pars(params, target=target)
    ## now get initial conds
    state_2 <- get_init(date0,date1,p2,int,slope,var)
    return(list(state=state_2,params=p2))
}    


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




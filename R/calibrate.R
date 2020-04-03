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
badness <- function(delta, params, target, state) {
    p_new <- adjust_params(delta[1],delta[2],params)
    s_new <- rep(NA,length(target))
    setNames(s_new,names(target))
    GIm <- get_GImoments(p)
    for (i in names(target)) {
        s_new[[i]] <- switch(i,
                              r=get_r(p_new),
                              R0=get_R0(p_new),
                              Gbar=GIm[["Gbar"]],
                              kappa=GIm[["kappa"]])
    }
    return(sum((s_new-target)^2))
}

## round-trip: should be zero
## badness(c(0,0), params=p1,target=c(r=get_r(p1,state), R0=get_R0(p1)),
##        state=state)

## 
## badness(c(0,0), params=p1,target=c(r=0.23,R0=3), state=state)

##' adjust pars to match targets
##' @export
fix_pars <- function(params, state, target=c(r=0.23,R0=3)) {
    opt1 <- optim(par=c(0,0), fn= badness, method="Nelder-Mead",
                  target=target, state=state, params=params)
    p_new <- adjust_params(opt1$par[1],opt1$par[2], params)
    return(p_new)
}

##' given a log-linear intercept and slope and target values for R0 or GI,
##' return initial conditions and parameters that will (approximately?)
##' recapitulate the observed log-linear dynamics
##' @param int log-linear intercept (expected value on date1)
##' @param slope log-linear slope (growth rate)
##' @param R0 desired/calibrated R0
##' @param GI desired/calibrated mean generation interval (stub)
##' @param date0 date for initial conditions
##' @param date1 initial date for regression
##' @param var variable for regression (H, D, ICU)
##' @export
calibrate <- function(int,slope,pop,params,R0=3,GI=NULL,
                            date0=ldmy("1-Mar-2020"),
                            date1=ldmy("25-Mar-2020"),
                            var="H") {
    ldmy <- lubridate::dmy  ## sugar) {
    ## first adjust base params to set r and R0 (or GI) to specified value
    params[["N"]] <- pop
    state <- make_state(N=pop,E0=1)
    p2 <- fix_pars(params,state,target=c(r=slope,R0=R0))
    ## now get initial conds
    list(state=get_init(date0,date1,p2,int,slope,var),
         params=p2)
}    


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
    s_new <- c(r=get_r(p_new,state),R0=get_R0(p_new))
    return((s_new[["r"]]-target[["r"]])^2 + (s_new[["R0"]] -target[["R0"]])^2)
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

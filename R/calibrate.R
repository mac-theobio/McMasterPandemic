
adjust_params <- function(log_delta,
                          pars_adj=list("beta0",
                                        c("sigma","gamma_s","gamma_m","gamma_a")),
                          params) {
    for (i in seq_along(log_delta)) {
        params[pars_adj[[i]]] <- params[pars_adj[[i]]]*exp(log_delta[i])
    }
    return(params)
}

## 
uniroot_target <- function(log_delta, params, target, pars_adj) {
    p_new <- adjust_params(log_delta, pars_adj=pars_adj, params)
    val <- switch(names(target),
                  r=get_r(p_new),
                  R0=get_R0(p_new),
                  Gbar=get_Gbar(p_new))
    return(val-target)
}


badness <- function(delta, params, target, pars_adj) {
    p_new <- adjust_params(delta, params, pars_adj=pars_adj)
    s_new <- rep(NA,length(target))
    names(s_new) <- names(target)
    GIm <- get_GI_moments(p_new)
    for (i in names(target)) {
        s_new[[i]] <- switch(i,
                             r=get_r(p_new),
                             R0=get_R0(p_new),
                             Gbar=get_Gbar(p_new),
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
##' @importFrom stats optim
##' @param params a parameter vector
##' @param target target values for one or more epidemic moments
##' @param pars_adj list of sets of parameters to adjust
##' @param u_interval interval for uniroot adjustment
##' @param r_method method for fixing r (brute-force exponential simulation or JD's kernel-based approach?)
##' @export
## FIXME: automatically choose default pars_adj on the basis of target?
##  better checking of pars_adj (should be a list of named vectors with
##  length equal to length of target
fix_pars <- function(params, target=c(r=0.23,Gbar=6),
                     pars_adj=list("beta0",
                                   c("sigma","gamma_s","gamma_m","gamma_a")),
                     r_method=c("expsim","rmult"),
                     u_interval = c(-3,3))
{
    r_method <- match.arg(r_method)
    if (identical(names(target),"r") && r_method=="rmult") {
        k <- transKernel(params)$foi
        rm <- rmult(k,target[["r"]])
        p_new <- params
        p_new[["beta0"]] <- p_new[["beta0"]]*rm
    } else {
        if (length(target)==1) {
            u <- uniroot(f=uniroot_target,interval=u_interval,
                         params=params, target=target, pars_adj=pars_adj)
            p_new <- adjust_params(u$root, params, pars_adj = pars_adj)
        } else {
            argList <- list(par=c(0,0), fn= badness, target=target, params=params, pars_adj=pars_adj)
            argList$method <- "Nelder-Mead"
            
            opt1 <- do.call(optim, argList)
            p_new <- adjust_params(opt1$par, params, pars_adj=pars_adj)
        }
    }
    return(p_new)
}

##' given a log-linear intercept and slope and target values for R0 or GI,
##' return initial conditions and parameters that will (approximately?)
##' recapitulate the observed log-linear dynamics
##' @inheritParams get_init
##' @param int log-linear intercept (expected value on date1)
##' @param slope log-linear slope (growth rate)
##' @param target vector of target statistics (allowed: R0, Gbar (mean generation interval), or kappa (CV^2 of generation interval)
##' @param pop population size
##' @param date0 date for initial conditions
##' @param date1 initial date for regression
##' @param var variable for regression (H, D, ICU)
##' @importFrom utils tail
##' @importFrom dplyr select pull
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
                      init_target=NULL,
                      date0="1-Mar-2020",
                      date1="25-Mar-2020",
                      var="H",
                      sim_args=NULL) {
    ## FIXME: refactor/rejigger more consistently
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
    if (is.null(init_target)) {
        state_2 <- get_init(date0,date1,p2,int,slope,var)
        ## IC is fully determined by E0 in this case
        p2[["E0"]] <- state_2[["E"]]
    } else {
        ## here (if it ever works) IC is more complex, includes
        ##  many non-zero compartments
        state_2 <- get_init(date0,date1,p2,var="H",init_target=init_target,
                            sim_args=sim_args)
    }
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
##' @param init_target value of focal variable to hit on date1
##' @param sim_args additional parameters to pass to \code{run_sim} if shooting
##' @examples
##' params <- read_params(system.file("params","ICU1.csv",package="McMasterPandemic"))
##' ## eigenvector projection
##' get_init(var="H",int=2.85,slope=0.176, param=params)
##' ## shooting
##' get_init(var="H",init_target=50, param=params)
##' @export
get_init <- function(date0=ldmy("1-Mar-2020"),
                     date1=ldmy("25-Mar-2020"),
                     params,
                     int=NULL,
                     slope=NULL,
                     var="H",
                     init_target=NULL,
                     sim_args=NULL) {
    ## analytical, back-projecting slope and intercept
    if (is.null(init_target)) {
        ## find eigenvectors of Jacobian (at disease-free eq)
        dom_vec <- get_evec(params)
        ## back-project modeled var from date1 to date0
        ## (predicted (modeled var) on date0
        int0 <- exp(int - as.numeric(date1-date0)*slope)
        ## scale state to expected (modeled var)
        state <- make_state(N=params[["N"]],E0=1) ## E0 will be replaced
        expected <- round(int0*dom_vec/dom_vec[[var]])
        expected <- expected[expected>0]
        state[names(expected)] <- expected
        ## not important, but adjust to keep same total pop
        state[["S"]] <- params[["N"]] - sum(expected)
    } else {
        ## brute force
        s0 <- make_state(params=params)
        ## assume we are just changing E0
        sim_args <- c(nlist(params,
                           state=make_state(params=params),
                           start_date=date0,
                           end_date=date1),
                     sim_args)
        ufun <- function(E0) {
            sim_args$state[["E"]] <- E0
            ## run simulation, aggregate
            r <- aggregate(do.call(run_sim,sim_args))
            return(tail(r[[var]],1)-init_target)
        }
        params[["E0"]] <- uniroot(ufun,interval=c(0.001,2000))$root
        state <- make_state(params=params)
        ## save("date0","date1","params","state",file="tmp.RData")
    }
    return(state)
}

##' estimate magnitude of drop in beta0 at known dates
##' @inheritParams calibrate
##' @inheritParams run_sim_break
##' @param optim_args additional args to \code{optim}
##' @param end_date ending date of sim
##' @param data data to match
##' @importFrom stats plogis
##' @export
## FIXME: take end_date from data
get_break <- function(date0=ldmy("1-Mar-2020"),
                      end_date=ldmy("8-Apr-2020"),
                      break_dates=c("20-Mar-2020"),
                      params,
                      data,
                      var="H",
                      sim_args=NULL,
                      ## FIXME, unused
                      optim_args=list())
{
    value <- NULL ## global var check
    ## FIXME: check order of dates?
    ## brute force
    s0 <- make_state(params=params)
    ## we are changing beta0 at a pre-specified set of break dates
    ## (by an unknown amount)
    mle_fun <- function(p,data) {
        ## inverse-link parameters
        rel_beta0 <- plogis(p[seq(nbrk)])
        nb_disp <- exp(p[-(1:seq(nbrk))])
        r <- do.call(run_sim_break,
                     c(nlist(params,
                             break_dates,
                             start_date=date0,
                             end_date,
                             rel_beta0),
                       sim_args))
        ## run simulation, aggregate
        r <- (aggregate(r, pivot=TRUE)
            %>% dplyr::rename(pred="value")
        )
        ## match up sim results with specified data
        names(data) <- tolower(names(data)) ## ugh
        data <- dplyr::filter(data,var %in% unique(r$var))
        r2 <- (dplyr::left_join(data,r,by=c("date","var"))
            %>% tidyr::drop_na(value,pred))
        ## FIXME: why do we have an NA in pred??
        ## compute negative log-likelihood
        ## FIXME assuming a single nb_disp for now
        ret <- with(r2,-sum(dnbinom(value,mu=pred,size=nb_disp,log=TRUE)))
        return(ret)
    }
    ## use optim to start with; maybe switch to mle2 later
    nbrk <- length(break_dates)
    nvar <- length(var)
    optim(par=rep(0,nbrk+nvar), fn=mle_fun, data=data)
}

##' run simulation with one or more breakpoints
## FIXME: make rel_beta0 part of params??? probably
## FIXME: roll into run_sim?
##' @param ... additional arguments to \code{run_sim}
##' @param params parameters
##' @param break_dates dates of breakpoints in transmission
##' @param rel_beta0 numeric vector (same length as \code{break_dates}): transmission relative to original value after each breakpoint
##' @export
run_sim_break <- function(params,
                          break_dates,
                          rel_beta0,
                          ...) {
    ## construct time-varying frame, parameters
    timevar <- data.frame(Date=break_dates,
                          Symbol="beta0",
                          Relative_value=rel_beta0)
    sim_args <- c(list(...),
                  nlist(params,
                        params_timevar=timevar,
                        state=make_state(params=params))
                  )
    do.call(run_sim,sim_args)
}
                        

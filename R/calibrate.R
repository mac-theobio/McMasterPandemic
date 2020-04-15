
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
##' @param debug debug?
##' @export
## FIXME: automatically choose default pars_adj on the basis of target?
##  better checking of pars_adj (should be a list of named vectors with
##  length equal to length of target, should be appropriate to target ...)
fix_pars <- function(params, target=c(r=0.23,Gbar=6),
                     pars_adj=list("beta0",
                                   c("sigma","gamma_s","gamma_m","gamma_a")),
                     r_method=c("expsim","rmult"),
                     u_interval = c(-3,3),
                     debug=FALSE)
{
    r_method <- match.arg(r_method)
    if (identical(names(target),"r") && r_method=="rmult") {
        k <- transKernel(params)$foi
        rm <- rmult(k,target[["r"]])
        p_new <- params
        p_new[["beta0"]] <- p_new[["beta0"]]*rm
    } else {
        if (length(target)==1) {
            if (debug) {
                cat("plotting uniroot curve\n")
                pvec <- seq(u_interval[1], u_interval[2], length.out=101)
                uvec <- vapply(pvec, uniroot_target, params=params,target=target, pars_adj=pars_adj, FUN.VALUE=numeric(1))
                plot(pvec,uvec)
                graphics::abline(v=0,lty=2)
            }
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
##' params <- read_params("ICU1.csv")
##' cc1 <- calibrate_slopeint(int=log(200),slope=0.23,pop=1e6,params)
##' summary(cc1$params)
##' cc2 <- calibrate_slopeint(int=log(200),slope=0.23,pop=1e6,params,target=c(R0=3))
##' summary(cc2$params)
## FIXME: warn if starting conditions are too low 
##' 
##' @export
calibrate_slopeint <- function(int,slope,pop,params,
                      target=c(Gbar=6),
                      init_target=NULL,
                      date0="1-Mar-2020",
                      date1="25-Mar-2020",
                      var="H",
                      sim_args=NULL) {
    ## FIXME: refactor/rejigger more consistently
    date0 <- anydate(date0); date1 <- anydate(date1)
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
get_init <- function(date0=anydate("1-Mar-2020"),
                     date1=anydate("25-Mar-2020"),
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
            r <- condense(do.call(run_sim,sim_args))
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
##' @param start_date starting date
##' @param end_date ending date of sim
##' @param data data to match
##' @param var variable to match
##' @importFrom stats plogis
##' @export
## FIXME: take end_date from data
get_break <- function(start_date=anydate("1-Mar-2020"),
                      end_date=anydate("8-Apr-2020"),
                      break_dates=c("20-Mar-2020"),
                      params,
                      data,
                      var="H",
                      sim_args=NULL,
                      ## FIXME, unused
                      optim_args=list())
{
    value <- pred <- NULL ## global var check
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
                             start_date=start_date,
                             end_date,
                             rel_beta0),
                       sim_args))
        ## run simulation, aggregate
        r <- (pivot(condense(r))
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
                        

##' Calibrate parameters 
##' @param base_params base parameters
##' @param target calibration values: mean generation interval, initial growth rate, initial value of calibration variable (\code{X0})
##' @param start_date starting date for simulations (should be far enough before the start of the data to allow sorting-out of initial conditions)
##' @param reg_date first day of calibration data used in log-linear regression
##' @param init_var variable to use in initial-condition calibration matching
##' @param sim_args list of additional arguments to pass to \code{\link{run_sim}}
##' @export
calib2 <- function(base_params,
                   target=c(Gbar=6,r=0.23,X0=200),
                   start_date,
                   reg_date,
                   init_var="H",
                   sim_args=NULL) {
    start_date <- anydate(start_date)
    reg_date <- anydate(reg_date)
    p <- fix_pars(base_params, target=c(Gbar=target[["Gbar"]]),
                  pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a")))
    stopifnot(all.equal(get_Gbar(p),target[["Gbar"]],tolerance=1e-4))
    ## fix params to *initial* slope (for now)
    p <- fix_pars(p, target=c(r=target[["r"]]),
                  pars_adj=list("beta0"),u_interval=c(-1,2))
    ## FIXME: problems with r-calibration here if using r_method="rmult"
    ## problem, or mismatch?
    stopifnot(all.equal(get_r(p),target[["r"]],tolerance=1e-4))
    ## calibrate initial state (shooting)
    state <- get_init(date0=start_date,date1=reg_date,p,var=init_var,
                      init_target=target[["X0"]], sim_args=sim_args)
    p[["E0"]] <- state[["E"]]
    return(p)
}

##' calibrate and run simulation
##' obsolete, replaced by forecast_sim?
##' @inheritParams calib2
##' @param end_date ending date for simulation
##' @param return_val return values: "sim"=full simulation (wide format); "aggsim"=aggregated and pivoted simulation output, "vals_only" = vector containing only values output from simulation
##' @export
sim_fun <- function(target, base_params,
                    start_date="1-Mar-2020",
                    reg_date,
                    end_date="1-May-2020",
                    init_var="H",
                    sim_args=NULL,
                    return_val=c("sim","aggsim","vals_only")) {
    return_val <- match.arg(return_val)
    ## calibration
    p <- calib2(base_params,
                target=target,
                start_date=start_date,
                reg_date=reg_date,
                init_var=init_var,
                sim_args=sim_args)
    state <- make_state(params=p)
    arg_list <- nlist(params=p,
                      state,
                      start_date,
                      end_date)
    r <- do.call(run_sim,arg_list)
    a <- pivot(condense(r))
    return(switch(return_val,
                  sim=r,
                  aggsim=a,
                  vals_only=a$value))
}

##' simulate based on a vector of parameters (including both time-varying change parameters, initial conditions, and other dynamical parameters), for fitting or forecasting
##' @importFrom stats update
##' @inheritParams calibrate
##' @param p vector of parameters
##' @param return_val specify values to return (aggregated simulation, or just the values?)
##' @export
forecast_sim <- function(p, opt_pars, base_params, start_date, end_date, break_dates,
                         fixed_pars = NULL,
                         sim_args=NULL, aggregate_args=NULL,
                         ## FIXME: return_val is redundant with sim_fun
                         return_val=c("aggsim","vals_only"))
{
    return_val <- match.arg(return_val)
    ## restructure and inverse-link parameters
    pp <- invlink_trans(restore(p, opt_pars, fixed_pars))
    ## substitute into parameters
    params <- update(base_params, params=pp$params)
    ## run simulation (uses params to set initial values)
    r <- do.call(run_sim_break,
                 c(nlist(params,
                         start_date,
                         end_date,
                         break_dates,
                         rel_beta0=pp$rel_beta0),
                   sim_args))
    ## aggregate
    r_agg <- condense(r)
    if (!is.null(aggregate_args)) {
        r_agg <- do.call(aggregate, c(list(r_agg),aggregate_args))
    }
    r_agg <- pivot(r_agg)
    ret <- switch(return_val,
                aggsim=r_agg,
                vals_only=r_agg$value
                )
    return(ret)
}

##' calibrate via negative binomial MLE, simultaneously fitting initial conditions, initial growth rate, time-changes in growth rate, and dispersion parameters
##' @param start_date starting date for sims (far enough back to allow states to sort themselves out)
##' @param start_date_offset days to go back before first data value
##' @param end_date ending date
##' @param break_dates specified breakpoints in beta0
##' @param base_params baseline parameters
##' @param data a data set to compare to, containing date/var/value (current version assumes that only a single state var is included)
##' @param opt_pars starting parameters (and structure).  Parameters that are part of the \code{params_pansim} parameter vector can be specified within the \code{params} element (with prefixes if they are transformed); other parameters can include distributional parameters or time-varying parameters
##' @param fixed_pars parameters to fix
##' @param sim_args additional arguments to pass to \code{\link{run_sim}}
##' @param aggregate_args arguments passed to \code{\link{aggregate.pansim}}
##' @param optim_args arguments passed to \code{\link{optim}}
##' @param debug print debugging messages?
##' @param debug_plot plot debugging curves?
##' @export
calibrate <- function(start_date=min(data$date)-start_date_offset,
                          start_date_offset=15,
                          end_date=max(data$date),
                          break_dates=c("23-Mar-2020","30-Mar-2020"),
                          base_params,
                          data,
                          opt_pars=list(params=c(log_E0=4,
                                                 log_beta0=-1),
                                        log_rel_beta0=c(-1,-1),
                                        log_nb_disp=0),
                          fixed_pars=NULL,
                          sim_args=NULL,
                          aggregate_args=NULL,
                          optim_args=NULL,
                          debug=FALSE,
                          debug_plot=FALSE)
                          
{
    if (debug) {
        cat("start date: ", format(start_date),
            "; end_date: ", format(end_date), "\n")
    }
    ## FIXME: check that appropriate var names are present
    ## translate state variables names in data to expected internal names (e.g. newConfirmations -> reports)
    ## FIXME: should this be external?
    data$var <- trans_state_vars(data$var)
    value <- pred <- NULL ## global var check
    ## FIXME: check order of dates?
    ## brute force
    ## we are changing beta0 at a pre-specified set of break dates
    ## (by an unknown amount)
    ## FIXME: allow for aggregation of reports etc.
    ##  ?? build this as an optional argument to
    ##   aggregate.pansim
    mle_fun <- function(p,data,do_plot=FALSE) {
        var <- NULL 
        r <- (do.call(forecast_sim,
                     nlist(p, opt_pars, base_params, start_date, end_date, break_dates,
                           sim_args, aggregate_args))
            %>% dplyr::rename(pred="value")
        )
        ## match up sim results with specified data
        ## FIXME: do these steps (names fix, filter) outside ?
        names(data) <- tolower(names(data)) ## ugh
        ## discard unused state variables
        data2 <- dplyr::filter(data,var %in% unique(r$var))
        ## keep only dates present in data
        r2 <- (dplyr::left_join(data2,r,by=c("date","var"))
            %>% tidyr::drop_na(value,pred))         ## FIXME: why do we have an NA in pred??
        ## compute negative log-likelihood
        ## FIXME assuming a single nb_disp for now
        if (do_plot) {
            plot(value~date, data=r2, log="y")
            with(r2,lines(date,pred))
        }
        ## need this for NB parameter
        pp <- invlink_trans(restore(p, opt_pars, fixed_pars))
        ret <- with(r2,-sum(dnbinom(value,mu=pred,size=pp$nb_disp,log=TRUE)))
        ## FIXME: add evaluation number?
        if (debug) cat(unlist(pp),ret,"\n")
        return(ret)
    }
    opt_inputs <- unlist(opt_pars)
    if (!is.null(fixed_pars)) {
        opt_inputs <- opt_inputs[setdiff(names(opt_inputs), names(unlist(fixed_pars)))]
    }
    ## use optim to start with; maybe switch to mle2 later
    opt <- do.call(optim,
            c(list(par=opt_inputs, fn=mle_fun, data=data,
                   do_plot=debug_plot),
              optim_args))
    if (opt$convergence>0) warning("convergence problem in optim()")
    attr(opt,"forecast_args") <- nlist(start_date,
                                       end_date,
                                       break_dates,
                                       base_params,
                                       opt_pars,
                                       fixed_pars,
                                       sim_args,
                                       aggregate_args)
    class(opt) <- c("fit_pansim",class(opt))
    return(opt)
}


##' find confidence envelopes by simulation
##' @param fit output from \code{calibrate}
##' @param nsim number of simulations
##' @param seed random-number seed
##' @param forecast_args arguments to pass to \code{forecast_sim}
##' @param qvec vector of quantiles
##' @param qnames quantile names
##' @export
## FIXME: way to add args to forecast_args list, e.g. stochastic components?
## FIXME: use bbmle::pop_pred_samp?
forecast_ensemble <- function(fit,
                              nsim=200,
                              forecast_args=attr(fit,"forecast_args"),
                              qvec=c(0.05,0.5,0.95),
                              qnames=c("lwr","value","upr"),
                              seed=NULL
                              ) {

    var <- NULL
    ## FIXME: include baseline forecast as well?
    
    ## parameter ensemble
    if (!is.null(seed)) set.seed(seed)
        ## HACK: for aggregated data, NB fit is unhappy because there is severe underdispersion (because of
        ##  overfitting to time series); NB disp parameter is >>> 1
        ##  we seem to be able to get away with ignoring it completely here
        ##  (not needed for forecast ...)
    ## might help to fix it ...
    ff <- function(p, return_val="vals_only") {
        do.call(forecast_sim,
                c(list(p,return_val=return_val),forecast_args))
    }

    ## baseline fit
    r <- ff(fit$par, return_val="aggsim")

    ## Wald sample
    ## FIXME: count number of distribution params
    parnum <- length(fit$par) - 1
    e_pars <- as.data.frame(MASS::mvrnorm(nsim,
                                          mu=fit$par[1:parnum],
                                          Sigma=solve(fit$hessian[1:parnum,1:parnum])))

    ## run for all param vals in ensemble
    ## tried with purrr::pmap but too much of a headache
    t1 <- system.time(e_res <- plyr::alply(as.matrix(e_pars)
                                         , .margins=1
                                         , .fun=ff))
    ## get quantiles per observation
    e_res2 <- (e_res %>% dplyr::bind_cols()
        %>% apply(1,stats::quantile,qvec,na.rm=TRUE)
        %>% t()
        %>% dplyr::as_tibble()
        %>% setNames(qnames)
    )
    
    ## date/var values
    e0 <- (dplyr::select(r,date,var)
        %>% dplyr::as_tibble()
    )
    ## combine quantiles with the original date/var columns
    e_res3 <- dplyr::bind_cols(e0, e_res2)
    return(e_res3)
}

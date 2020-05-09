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


## FIXME: originally called sim_fun, maybe used by some obsolete code
## its existence as sim_fun() confuses byte compiler ...
##' calibrate and run simulation
##' obsolete, replaced by forecast_sim?
##' @inheritParams calib2
##' @param end_date ending date for simulation
##' @param return_val return values: "sim"=full simulation (wide format); "aggsim"=aggregated and pivoted simulation output, "vals_only" = vector containing only values output from simulation
##' @export
sim_fun_old <- function(target, base_params,
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

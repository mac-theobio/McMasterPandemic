##' simulate based on a vector of parameters (including both time-varying change parameters, initial conditions, and other dynamical parameters), for fitting or forecasting
##' @inheritParams get_break_gen
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
    params <- update(base_params, E0=pp[["E0"]], beta0=pp[["beta0"]])
    ## run simulation (uses params to set initial values)
    r <- do.call(run_sim_break,
                 c(nlist(params,
                         start_date,
                         end_date,
                         break_dates,
                         rel_beta0=pp$rel_beta0),
                   sim_args))
    ## aggregate
    r_agg <- do.call(aggregate,
                     c(list(r, pivot=TRUE),
                       aggregate_args))
    ret <- switch(return_val,
                aggsim=r_agg,
                vals_only=r_agg$value
                )
    return(ret)
}

##' calibrate via negative binomial MLE, simultaneously fitting initial conditions, initial growth rate, time-changes in growth rate, and dispersion parameters
##' @param start_date starting date for sims (far enough back to allow states to sort themselves out)
##' @param end_date ending date
##' @param break_dates specified breakpoints in beta0
##' @param base_params baseline parameters
##' @param data a data set to compare to, containing date/var/value (current version assumes that only a single state var is included)
##' @param var variable to compare to (FIXME: this should be inferred from data!
##' @param opt_pars starting parameters (and structure)
##' @param sim_args additional arguments to pass to \code{\link{run_sim}}
##' @param aggregate_args arguments passed to \code{\link{aggregate.pansim}}
##' @param optim_args arguments passed to \code{\link{optim}}
##' @export
get_break_gen <- function(start_date=min(data$date)-start_date_offset,
                          start_date_offset=15,
                          end_date=max(data$date),
                          break_dates=c("23-Mar-2020","30-Mar-2020"),
                          base_params,
                          data,
                          var="report",
                          opt_pars=list(log_E0=4,
                                        log_beta0=-1,
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
    do.call(optim,
            c(list(par=opt_inputs, fn=mle_fun, data=data,
                   do_plot=debug_plot),
              optim_args))
}


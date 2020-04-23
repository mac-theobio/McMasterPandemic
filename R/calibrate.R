
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

##' run simulation with one or more breakpoints
## FIXME: make rel_beta0 part of params??? probably
## FIXME: roll into run_sim?
##' @param ... additional arguments to \code{run_sim}
##' @param params parameters
##' @param break_dates dates of breakpoints in transmission
##' @param rel_beta0 numeric vector (same length as \code{break_dates}): transmission relative to original value after each breakpoint
##' @examples
##' params <- read_params("ICU1.csv")
##' r1 <- run_sim_break(params, break_dates="2020-03-01",
##'                    start_date="2020-02-01", end_date="2020-04-01",
##'                    rel_beta0 = 0.5)
##' plot(r1,log=TRUE)
##' r2 <- run_sim_break(params, 
##'                    start_date="2020-02-01", end_date="2020-04-01")
##' plot(r2,log=TRUE)
##' @export
run_sim_break <- function(params,
                          break_dates=NULL,
                          rel_beta0,
                          ...) {
    sim_args <- c(list(...),
                  nlist(params,
                        state=make_state(params=params)))
    if (!is.null(break_dates)) {
        ## construct time-varying frame, parameters
        timevar <- data.frame(Date=break_dates,
                              Symbol="beta0",
                              Relative_value=rel_beta0)
        sim_args <- c(sim_args,
                      list(params_timevar=timevar))
    }
    do.call(run_sim,sim_args)
}

##' simulate based on a vector of parameters (including both time-varying change parameters, initial conditions, and other dynamical parameters), for fitting or forecasting
##' @importFrom stats update
##' @inheritParams calibrate
##' @param p vector of parameters
##' @param condense_args arguments to pass to \code{\link{condense}} (via \code{\link{run_sim}})
##' @param stoch stochastic settings (see \code{\link{run_sim}})
##' @param return_val specify values to return (aggregated simulation, or just the values?)
##' @export
forecast_sim <- function(p, opt_pars, base_params, start_date, end_date, break_dates,
                         fixed_pars = NULL,
                         stoch = NULL,
                         sim_args=NULL, aggregate_args=NULL, condense_args=NULL,
                         ## FIXME: return_val is redundant with sim_fun
                         return_val=c("aggsim","vals_only"),
                         debug = FALSE)
{
    return_val <- match.arg(return_val)
    if (!is.null(stoch)) {
        if (is.null(sim_args)) {
            sim_args <- list(stoch=stoch)
        } else {
            sim_args <- c(sim_args, list(stoch=stoch))
        }
    }
    ## restructure and inverse-link parameters
    pp <- invlink_trans(restore(p, opt_pars, fixed_pars))
    ## substitute into parameters
    params <- update(base_params, params=pp$params, .list=TRUE)
    if (debug) cat("forecast ",params[["beta0"]],"\n")
    ## run simulation (uses params to set initial values)
    r <- do.call(run_sim_break,
                 c(nlist(params,
                         start_date,
                         end_date,
                         break_dates,
                         rel_beta0=pp$rel_beta0,
                         condense_args),
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

mle_fun <- function(p, data, debug=FALSE, debug_plot=FALSE,
                    opt_pars, base_params, start_date, end_date, break_dates=NULL,
                    sim_args=NULL, aggregate_args=NULL,
                    priors=NULL, ...) {
    ## ... is to drop any extra crap that gets in there
    ## opt_pars <- base_params <- start_date <- end_date <- NULL
    ## break_dates <- sim_args <- aggregate_args <- NULL
    if (debug) cat(p,"\n")
    var <- pred <- value <- NULL 
    r <- (do.call(forecast_sim,
                  nlist(p, opt_pars, base_params, start_date, end_date, break_dates,
                        sim_args, aggregate_args, debug))
        %>% dplyr::rename(pred="value")
    )
    ## ggplot(r,aes(date,pred,colour=var)) + geom_line() + scale_y_log10() + geom_point(data=data,aes(y=value))
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
    if (debug_plot) {
        vv <- unique(r2$var)
        if (length(vv)==1) {
            suppressWarnings(plot(value~date, data=r2, log="y"))
            with(r2,lines(date,pred))
        } else {
            with(r2,plot(date,value,col=as.numeric(factor(var)),log="y"))
            r2s <- split(r2,r2$var)
            invisible(Map(function(x,c) lines(x$date,x$pred, col=c), r2s, seq_along(r2s)))
        }
    }
    ## need this for NB parameter
    ## FIXME: fixed params can now be handled through mle2?
    pp <- invlink_trans(restore(p, opt_pars))
    dvals <- with(r2,dnbinom(value,mu=pred,size=pp$nb_disp,log=TRUE))
    ## clamp NaN/NA values to worst obs
    dvals[is.na(dvals)] <- min(dvals,na.rm=TRUE)
    ret <- -sum(dvals)
    if (!is.null(priors)) {
        browser()
    }
    ## FIXME: add evaluation number?
    if (debug) cat(unlist(pp),ret,"\n")
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
##' @param mle2_control control args for mle2
##' @param mle2_method method arg for mle2
##' @param mle2_args additional arguments for mle2
##' @param debug print debugging messages?
##' @param debug_plot plot debugging curves?
##' @importFrom graphics lines
##' @importFrom bbmle parnames<- mle2
## DON'T import stats::coef !
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
                      debug=FALSE,
                      debug_plot=FALSE,
                      mle2_method="Nelder-Mead",
                      mle2_control=list(maxit=10000),
                      mle2_args=list())
{
    cc <- match.call()
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
    ## FIXME: allow for aggregation of reports etc.
    ## MLE
    opt_inputs <- unlist(opt_pars)
    ## n.b.: this has to be hacked dynamically ...
    parnames(mle_fun) <- names(opt_inputs)
    mle_data <- nlist(start_date,
          end_date,
          break_dates,
          base_params,
          opt_pars,
          sim_args,
          fixed_pars,
          aggregate_args,
          debug,
          debug_plot,
          data)
    ## unnecessary
    ## if (utils::packageVersion("bbmle")<"1.0.23.2") stop("please remotes::install('bbolker/bbmle') to get newest version")
    opt <- do.call(bbmle::mle2,
                   c(list(minuslogl=mle_fun
                        , start=opt_inputs
                        , data=mle_data
                        , vecpar=TRUE
                        , method=mle2_method
                        , control=mle2_control
                        ## , namedrop_hack = FALSE  ## unnecessary
                          ),
                     mle2_args))
    res <- list(mle2=opt,
                forecast_args=nlist(start_date,
                                    end_date,
                                    break_dates,
                                    base_params,
                                    opt_pars,
                                    fixed_pars,
                                    sim_args,
                                    aggregate_args),
                call=cc)
    class(res) <- "fit_pansim"
    return(res)
}

##' find confidence envelopes by simulation
##' @param fit output from \code{calibrate}
##' @param nsim number of simulations
##' @param seed random-number seed
##' @param forecast_args arguments to pass to \code{forecast_sim}
##' @param qvec vector of quantiles
##' @param qnames quantile names
##' @param .progress progress bar?
##' @export
## FIXME: way to add args to forecast_args list, e.g. stochastic components?
## FIXME: use bbmle::pop_pred_samp?
forecast_ensemble <- function(fit,
                              nsim=200,
                              forecast_args=fit$forecast_args,
                              qvec=c(0.05,0.5,0.95),
                              qnames=c("lwr","value","upr"),
                              seed=NULL,
                              .progress=if (interactive()) "text" else "none"
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
    r <- ff(coef(fit$mle2), return_val="aggsim")

    ## Wald sample
    ## FIXME: count number of distribution params?
    ## FIXME: use pop_pred_samp()?
    e_pars <- as.data.frame(MASS::mvrnorm(nsim,
                                          mu=coef(fit$mle2),
                                          Sigma=bbmle::vcov(fit$mle2)))

    ## run for all param vals in ensemble
    ## tried with purrr::pmap but too much of a headache
    t1 <- system.time(e_res <- plyr::alply(as.matrix(e_pars)
                                         , .margins=1
                                         , .fun=ff
                                         , .progress=.progress  ))
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

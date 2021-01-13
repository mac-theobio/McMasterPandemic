
## utility for fix_pars: adjust specified set of parameters by a specified amount
adjust_params <- function(log_delta,
                          pars_adj=list("beta0",
                                        c("sigma","gamma_s","gamma_m","gamma_a")),
                          params) {
    for (i in seq_along(log_delta)) {
        params[pars_adj[[i]]] <- params[pars_adj[[i]]]*exp(log_delta[i])
    }
    return(params)
}

## utility for fix_pars: define an objective function for uniroot()
uniroot_target <- function(log_delta, params, target, pars_adj) {
    p_new <- adjust_params(log_delta, pars_adj=pars_adj, params)
    val <- switch(names(target),
                  r=get_r(p_new),
                  R0=get_R0(p_new),
                  Gbar=get_Gbar(p_new))
    return(val-target)
}


## objective function for fix_pars if we are trying to adjust multiple parameters simultaneously
## (i.e. we minimize sum of squares of the *vector* of deviations from target for all the parameters)
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
##' @examples
##' pp <- read_params("ICU1.csv")
##' summary(pp)
##' pp2 <- fix_pars(pp,debug=TRUE)
##' summary(pp2)
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
            if (debug) {
                argList$control  <- list(trace=1)
            }
            ## with(argList,fn(par,params=params,pars_adj=pars_adj,target=target))
            opt1 <- do.call(optim, argList)
            p_new <- adjust_params(opt1$par, params, pars_adj=pars_adj)
        }
    }
    return(p_new)
}


##' run with log-linear model applied to beta
##' @examples
##' params <- read_params("ICU1.csv")
##' ## UNFINISHED
##' @inheritParams run_sim_break
##' @export
run_sim_loglin <- function(params,
                           extra_pars=NULL,
                           time_args=NULL,
                           sim_args=list(),
                           return_timevar=FALSE,
                           ...) {
    ## FIXME: DRY from run_sim_break/run_sim_mobility
    X <- X_date <- NULL  ## variable checker
    unpack(time_args)
    ## X: model matrix
    ## X_date: dates corresponding to rows of model matrix
    ## 
    ## strip 
    time_beta <- extra_pars$time_beta
    extra_pars$time_beta <- NULL
    timevar <- NULL
    if (length(time_beta)>0) {  ## non-trivial model
        ## construct time-varying frame, parameters
        ## FIXME:: get rid of anydate() ?
        timevar <- dfs(Date=anydate(X_date),
                       Symbol="beta0",
                       Relative_value=exp(X %*% time_beta))  ## log-linear model for beta
    }

    if (!is.null(timevar) && length(sim_args)>0) {
        ## FIXME: under what circumstances can sim_args be empty?
        ## should we check start_date and end_date separately?
        if (sim_args$start_date < min(timevar$Date)) {
            constant_rel_val <- data.frame(Date= seq.Date(as.Date(sim_args$start_date),as.Date(min(timevar$Date))-1, by =1)
                                           , Symbol = "beta0"
                                           , Relative_value = 1)
            timevar <- rbind(constant_rel_val,timevar)
        }
        if (sim_args$end_date > max(timevar$Date)){
            freeze_dat <- data.frame(Date = seq.Date(as.Date(max(timevar$Date))+1, as.Date(sim_args$end_date), by = 1)
                                     , Symbol = "beta0"
                                     , Relative_value = timevar$Relative_value[nrow(timevar)]
            )
            timevar <- rbind(timevar,freeze_dat)
        }
    }
    ## This is a temp fix. Need to look inside run_sim time switches at these discontinuities rather than stitching bt
    
    
    if ("testing_data" %in% names(time_args)) {
        freeze_testing <- data.frame()
        if (sim_args$end_date > max(time_args$testing_data$Date)){
            freeze_testing <- data.frame(Date = seq.Date(as.Date(max(time_args$testing_data$Date))+1, as.Date(sim_args$end_date), by = 1)
                                     , Symbol = "testing_intensity"
                                     , Relative_value = time_args$testing_data$Relative_value[nrow(time_args$testing_data)]
            )
        }
        testing_dat <- rbind(time_args$testing_data,freeze_testing)
        timevar <- rbind(timevar,testing_dat)
    }
    if (return_timevar) return(timevar)
    sim_args <- c(sim_args
                , extra_pars
                , nlist(params,
                        state=make_state(params=params),
                        params_timevar=timevar))
    do.call(run_sim,sim_args)
}

##' run with transmission propto relative mobility 
##' @inheritParams run_sim_break
##' @export
run_sim_mobility <- function(params,
                             extra_pars=NULL,
                             time_args=NULL,
                             sim_args=list(),
                             return_timevar=FALSE,
                             ...) {
    ## FIXME: DRY from run_sim_break
    mob_value <- mob_startdate <- NULL
    unpack(time_args)
    ## mob_value: numeric vector of relative mobility per day (starting at 1 and decreasing)
    ## mob_startdate: starting date for the relative mobility sequence (anything before this == 1 rel mobility)
    ## strip 
    mob_power <- extra_pars$mob_power
    extra_pars$mob_power <- NULL
    ## construct time-varying frame, parameters
    timevar <- dfs(Date=seq.Date(from=mob_startdate,
                                        length.out=length(mob_value),
                                        by="1 day"),
                          Symbol="beta0",
                   Relative_value=mob_value^mob_power)
    if (return_timevar) return(timevar)
    sim_args <- c(sim_args
                , extra_pars
                , nlist(params,
                        state=make_state(params=params),
                        params_timevar=timevar))
    do.call(run_sim,sim_args)
}
    
##' run simulation with one or more breakpoints
## FIXME: make rel_beta0 part of params???
## FIXME: roll into run_sim?
## FIXME: generalize
##' @param ... additional arguments to \code{run_sim}
##' @param params parameters
##' @param time_args list containing \code{break_dates}
##' @param sim_args parameters to pass to \code{\link{run_sim}}
##' @param extra_pars parameters that are used to set up time-varying parameters, etc., but \emph{not} used by \code{run_sim}
##' @param break_dates obsolete
##' @param return_timevar return data frame of beta by time?
##' @examples
##' params <- read_params("ICU1.csv")
##' r1 <- run_sim_break(params, time_args=list(break_dates="2020-03-01"),
##'                    sim_args=list(start_date="2020-02-01", end_date="2020-04-01"),
##'                    extra_pars=list(rel_beta0 = 0.5))
##' plot(r1,log=TRUE)
##' ## can also use it to run without breaks ...
##' r2 <- run_sim_break(params, sim_args=list(start_date="2020-02-01", end_date="2020-04-01"))
##' plot(r2,log=TRUE)
##' @export
run_sim_break <- function(params,
                          extra_pars=NULL,
                          time_args=NULL,
                          break_dates=NULL,
                          sim_args=list(),
                          return_timevar=FALSE,
                          ...) {
    if (!is.null(break_dates)) {
        warning("use of break_dates as a top-level parameter is deprecated: please use time_args=list(break_dates=...)")
        time_args <- list(break_dates=break_dates)
    }
    ## FIXME:: HACK for now
    ## other_args <- other_args[!grepl("nb_disp",names(other_args))]
    sim_args <- c(sim_args,
                  nlist(params,
                        state=make_state(params=params)))
    if (length(time_args)==1 && is.null(names(time_args))) {
        ## HACK:: namedrop() problems in mle2????
        names(time_args) <- "break_dates"
    }
    if (!is.null(time_args$break_dates)) {
        ## construct time-varying frame, parameters
        timevar <- data.frame(Date=anydate(time_args$break_dates),
                              Symbol="beta0",
                              Relative_value=extra_pars$rel_beta0)
        if (return_timevar) return(timevar)
        sim_args <- c(sim_args,
                      list(params_timevar=timevar))
    }
    do.call(run_sim,sim_args)
}

## keep rel_beta0 as time-varying params argument (now misnamed)
## PRIORS final = 0.2 - 0.5  plogis norm(mean=-0.75,sd=0.75/2)
## initial = 1
## date = day 63 (= 2020-03-19)
## scale: 1.9 to 2.3 days  ( Norm (2.1, sd=0.1))
run_sim_decay <- function(params,
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

##' simulate/forecast a single trajectory
##' 
##' simulate based on a vector of parameters (including both time-varying change parameters, initial conditions, and other dynamical parameters), for fitting or forecasting
##'
##' @importFrom stats update
##' @inheritParams calibrate
##' @inheritParams run_sim
##' @param p vector of parameters - on the link (log/logit) scale as appropriate; these are the parameters that will be adjusted during calibration
##' @param condense_args arguments to pass to \code{\link{condense}} (via \code{\link{run_sim}})
##' @param stoch stochastic settings (see \code{\link{run_sim}})
##' @param return_val specify values to return (aggregated simulation, or just the values?)
##' @param calc_Rt calculate and include R(t) in prediction/forecast?
##' @param ... extra args (ignored)
##' @examples
##' ff <- ont_cal1$forecast_args
##' op <- ff$opt_pars
##' p <- unlist(op)
##' params <- fix_pars(read_params("ICU1.csv"))
##' forecast_sim(p, op, base_params=params,ff$start_date, ff$end_date,
##'     time_args=ff$time_args)
##' @export
forecast_sim <- function(p, opt_pars, base_params, start_date, end_date,
                         time_args = NULL,
                         fixed_pars = NULL,
                         stoch = NULL,
                         stoch_start = NULL,
                         sim_args=list(),
                         aggregate_args=NULL,
                         ##condense = TRUE,
                         condense_args=NULL,  ## FIXME: ???
                         ## FIXME: return_val is redundant with sim_fun
                         return_val=c("aggsim","vals_only"),
                         sim_fun=run_sim_break,
                         calc_Rt = FALSE,
                         debug = FALSE,
                         params_timevar,
                         ...)

{
    ## FIXME: shouldn't need ..., but catches e.g. debug_hist passed
    ## through?
    S <- Symbol <- rel_beta0 <- hetS <- zeta <- NULL ## global var checking
    return_val <- match.arg(return_val)
    sim_args <- c(sim_args,nlist(start_date, end_date))
    if (!is.null(stoch)) {
        sim_args <- c(sim_args, nlist(stoch, stoch_start))
    }
    ## restructure and inverse-link parameters
    ## if (debug) cat("forecast 0",opt_pars[["log_beta0"]],"\n")
    pp <- invlink_trans(restore(p, opt_pars, fixed_pars))
    ## substitute into parameters
    params <- update(base_params, params=pp$params, .list=TRUE)
    ## now drop params and dispersion parameters from pp, send everything else to sim_fun()
    pp <- pp[!grepl("^params$|nb_disp",names(pp))]
    ## if (debug) cat("forecast ",params[["beta0"]],"\n")
    ## run simulation (uses params to set initial values)
    all_sim_args <- c(nlist(params,
                         extra_pars=pp,
                         start_date,
                         end_date,
                         time_args,
                         condense_args,
                         sim_args,
                         params_timevar),
                      pp)
    r <- do.call(sim_fun, all_sim_args)
    ## FIXME: remove? already condensed?
    ## if (condense) r_agg <- condense(r)
    r_agg <- r
    ## aggregate
    if (!is.null(aggregate_args)) {
        r_agg <- do.call(aggregate, c(list(r_agg),aggregate_args))
    }
    if (calc_Rt) {
        R0_base <- summary(params)[["R0"]]
        ## retrieve time-varying beta
        bb <- do.call(sim_fun,c(all_sim_args,list(return_timevar=TRUE)))
        if (!is.null(bb)) {
            bb <- (bb
                %>% as_tibble()
                %>% dplyr::filter(Symbol=="beta0")
                %>% dplyr::select(-Symbol)
                %>% rename(rel_beta0="Relative_value",date="Date")
            )
            vars <- c("date","S")
            if (has_zeta(params)) vars <- c(vars,"hetS")
            x2 <- (full_join(bb,select(r_agg,vars),
                             by="date")
                %>% arrange(date))
        }  else {
            x2 <- r_agg %>% select(date,S, hetS) %>% mutate(rel_beta0=1)
        }
        x3 <- (x2
            %>% mutate_at("rel_beta0", fill_edge_values)
            %>% transmute(date=date, hetS=hetS, Rt=R0_base*rel_beta0)
        )
        if (has_zeta(params)) {
            x3 <- (x3
                %>% mutate_at("Rt", ~.*hetS)
                %>% select(-hetS)
            )
        }
        r_agg <- full_join(r_agg, x3, by="date")
    } ## calc_Rt
    r_agg <- tidyr::pivot_longer(r_agg,-date, names_to="var")
    ret <- switch(return_val,
                aggsim=r_agg,
                vals_only=r_agg$value
                )
    return(ret)
}

## utility: debugging plot
##' @importFrom graphics legend
do_debug_plot <- function(r2) {
    value <- NULL ## global var check
    vv <- unique(r2$var)
    if (length(vv)==1) {
        suppressWarnings(plot(value~date, data=r2, log="y", main=vv))
        with(r2,lines(date,pred))
    } else {
        with(subset(r2,value>0),
             plot(date,value,col=as.numeric(factor(var)),log="y"))
        r2s <- split(r2,r2$var)
        Map(function(x,c) lines(x$date,x$pred, col=c), r2s, seq_along(r2s))
        labs <- sort(unique(r2$var)) ## sort to match factor order
        legend("topright",col=seq(length(labs)),lty=1,legend=labs)
    }
}

##' negative log-likelihood function
##' @param p parameter vector (in unlisted form)
##' @param ... unused (but useful in case junk needs to be discarded)
##' @param checkpoint save file containing call information?
##' @param na_penalty value to add to NLL for NA values in log-likelihood
##' @inheritParams calibrate
##' @examples
##' library(dplyr)
##' p <- read_params("ICU1.csv")
##' op <- get_opt_pars(p)
##' dd <- ont_all %>% trans_state_vars() %>% filter(var %in% c("H","death"))
##' mle_fun(p=unlist(op), dd, opt_pars=op, base_params=p)
##' op <- op["params"] ## exclude log_nb_disp
##' try(mle_fun(p=unlist(op), dd, opt_pars=op, base_params=p))
##' p2 <- update(p, obs_disp=2)
##' mle_fun(p=unlist(op), dd, opt_pars=op, base_params=p2)
##' p3 <- update(p, obs_disp_H=2, obs_disp_death=2)
##' mle_fun(p=unlist(op), dd, opt_pars=op, base_params=p3)
##' @export
mle_fun <- function(p, data,
                    debug=FALSE,
                    debug_plot=FALSE,
                    debug_hist=FALSE,
                    opt_pars,
                    base_params,
                    start_date=min(data$date),
                    end_date=max(data$date),
                    time_args=NULL,
                    sim_args=NULL,
                    sim_fun=run_sim_break,
                    checkpoint=FALSE,
                    aggregate_args=NULL,
                    priors=NULL,
                    na_penalty=1e6,
                    ...) {

    ## browser()
    ## ... is to drop any extra crap that gets in there
    if (debug) cat("mle_fun: ",p,"\n")
    var <- pred <- value <- NULL    ## defeat global-variable checkers
    ## pass ALL of opt_pars and p to forecaster
    ##  nb_disp stuff isn't necessary within forecaster but we need p and opt_pars to stay in sync!
    f_args <- nlist(p
                  , opt_pars
                  , base_params
                  , start_date
                  , end_date
                  , time_args
                  , sim_args
                  , aggregate_args
                  , debug
                  , sim_fun)
    if (checkpoint) saveRDS(f_args,file=".mle_checkpoint.rds")
    r <- (do.call(forecast_sim, f_args)
        %>% dplyr::rename(pred="value") ## rename trajectory value to 'pred' (predicted)
    )
    ## ggplot(r,aes(date,pred,colour=var)) + geom_line() + scale_y_log10() + geom_point(data=data,aes(y=value))
    ## match up sim results with specified data
    ## FIXME: do these steps (names fix, filter) outside ?
    names(data) <- tolower(names(data)) ## ugh
    ## discard unused state variables
    data2 <- dplyr::filter(data,var %in% unique(r$var))
    ## join data and trajectory; keep only date/var combs present in data [NOT trajectory]
    r2 <- dplyr::left_join(data2,r,by=c("date","var")) %>% tidyr::drop_na(value)
    ## compute negative log-likelihood
    if (debug_plot) {
        do_debug_plot(r2)
    }
    ## need to reconstruct parameter set to get NB parameter(s)
    ## FIXME: fixed params can now be handled through mle2?
    pp <- invlink_trans(restore(p, opt_pars))
    nbd <- pp$nb_disp
    v <- unique(r2$var)
    if (is.null(nbd)) {
        if (!any(grepl("obs_disp",names(base_params)))) {
            stop("dispersion params must be specified in opt_pars *or* parameters")
        }
        ## dispersion parameters not specified in opt_pars ...
        if (length(vn <- grep("obs_disp_", names(base_params), value=TRUE))>0) {
            ## variable-specific dispersion specified
            nbd <- base_params[vn]
            names(nbd) <- gsub("obs_disp_","",vn)
        } else {
            nbd <- base_params[["obs_disp"]]
            ## will be replicated/named in next step
        }
    }
    if (is.null(names(nbd))) {
        if (length(nbd)>1) stop("nbdisp must be named or length==1")
        nbd <- setNames(rep(nbd,length(v)),v)
    }
    r2 <- merge(r2,data.frame(var=names(nbd),nb_disp=nbd), by="var")
    ## FIXME: more principled solution to dnbinom warnings etc.?
    ##  dnbinom wrapper that napredicts NAs?
    ## FIXED nb_disp hack, don't need to exponentiate any more ...
    dvals <- with(r2,-1*dnbinom(value,mu=pred,size=nb_disp,log=TRUE))
    ## clamp non-finite (Inf/NaN/NA) values to neg log likelihood of worst obs
    dvals[!is.finite(dvals)] <- max(dvals[is.finite(dvals)])+na_penalty
    ret <- sum(dvals)
    if (!is.null(priors)) {
        for (pr in priors) {
            pr <- pr[[2]] ## drop tilde
            pr <- add_d_log(pr)
            ret <- ret - eval(pr,list2env(pp))
        }
    }
    ## FIXME: add evaluation number?
    if (debug) {
        cat(unlist(pp),ret,"\n")
    }
    if (debug_hist) {
        update_debug_hist(unlist(pp),ret)
    }
    if (!is.finite(ret)) ret <- Inf  ## DEoptim hack
    return(ret)
}

debug_env <- new.env()
update_debug_hist <- function(params, NLL) {
    ## see https://stackoverflow.com/questions/63432138/equivalent-of-within-attach-etc-for-working-within-an-environment/63432196#63432196 for discussion
    debug_env$params <- params
    debug_env$NLL <- NLL
    with(debug_env,
    {
        if (debug_ctr>nrow(history_mat)) {
            history_mat <- rbind(history_mat, base_mat)  ## extend matrix
        }
        history_mat[debug_ctr,  ] <- c(params, NLL)
        debug_ctr <- debug_ctr + 1
    })
}

##' Estimate model parameters from data
##'
##' Given time series data and a set of starting
##' parameters/structure/etc., calibrate the model via negative
##' binomial maximum likelihood estimation (MLE), simultaneously
##' fitting initial conditions, initial growth rate, time-changes in
##' growth rate, and dispersion parameters.
##'
##' \code{\link[bbmle]{mle2}}
##' is used to estimate parameters by trajectory matching.
##' Differential evolution optimization is conducted
##' via \code{\link[DEoptim]{DEoptim}}.
##' 
##' @param start_date starting date for sims (far enough back to allow
##'     states to sort themselves out)
##' @param start_date_offset days to go back before first data value
##' @param end_date ending date
##' @param time_args list containing \code{break_dates} or other
##'     information needed for time-dependent variation
##' @param base_params baseline parameters (an object (vector?) of
##'     type \code{params_pansim} containing all of the parameters
##'     needed for a simulation; some may be overwritten during the
##'     calibration process)
##' @param data a data set to compare to, containing date/var/value
##'     (current version assumes that only a single state var is
##'     included)
##' @param opt_pars starting parameters (and structure).  Parameters
##'     that are part of the \code{params_pansim} parameter vector can
##'     be specified within the \code{params} element (with prefixes
##'     if they are transformed); other parameters can include
##'     distributional parameters or time-varying parameters
##' @param fixed_pars parameters to fix
##' @param sim_args additional arguments to pass to
##'     \code{\link{run_sim}}
##' @param aggregate_args arguments passed to
##'     \code{\link{aggregate.pansim}}
##' @param time_args arguments passed to \code{sim_fun}
##' @param break_dates legacy
##' @param mle2_control control args for mle2
##' @param mle2_method method arg for mle2
##' @param mle2_args additional arguments for mle2
##' @param debug print debugging messages?
##' @param debug_hist keep information on parameter history?
##' @param debug_plot plot debugging curves? (doesn't work with
##'     parallel DEoptim)
##' @param last_debug_plot plot debugging curve for \emph{only} last
##'     parameter set (stored in \code{.debug_plot.pdf} in current
##'     directory)
##' @param priors a list of tilde-delimited expressions giving prior
##'     distributions expressed in terms of the elements of
##'     \code{opt_pars},
##'     e.g. \code{list(~dlnorm(rel_beta0[1],meanlog=-1,sd=0.5))}
##' @param seed random-number seed (for DE)
##' @param use_DEoptim use differential evolution as first stage?
##' @param DE_args arguments for \code{\link{DEoptim}}
##' @param DE_lwr lower bounds for DE optimization
##' @param DE_upr upper bounds, ditto
##' @param DE_cores number of parallel workers for DE
##' @param condense_args arguments to pass to \code{\link{condense}}
##'     (via \code{\link{run_sim}}) [not implemented yet?]
##' @param sim_fun function for simulating a single run
##'     (e.g. \code{\link{run_sim_break}},
##'     \code{\link{run_sim_mobility}})
##' @importFrom graphics lines
##' @importFrom bbmle parnames<- mle2 
# DON'T import stats::coef !
##' @examples
##' library(dplyr)
##' params <- fix_pars(read_params("ICU1.csv"))
##'  opt_pars <- list(params=c(log_E0=4, log_beta0=-1,
##'         log_mu=log(params[["mu"]]), logit_phi1=qlogis(params[["phi1"]])),
##'                                   logit_rel_beta0=c(-1,-1),
##'                                    log_nb_disp=NULL)
##' dd <- (ont_all %>% trans_state_vars() %>% filter(var %in% c("report", "death", "H")))
##' \dontrun{
##'    cal1 <- calibrate(data=dd, base_params=params, opt_pars=opt_pars, debug_plot=TRUE)
##'   cal1_DE <- calibrate(data=dd, base_params=params, opt_pars=opt_pars, use_DEoptim=TRUE, DE_cores=1)
##'   cal2 <- calibrate(data=dd, base_params=params, opt_pars=opt_pars, use_DEoptim=TRUE)
##' 
##'    if (require(bbmle)) {
##'      -logLik(cal2$mle2)
##'    }
##'    attr(cal2,"de")$optim$bestval
##' }
##' @export
##' @importFrom stats var vcov
calibrate <- function(start_date=min(data$date)-start_date_offset,
                      start_date_offset=15,
                      end_date=max(data$date),
                      time_args=list(
                          break_dates=c("2020-Mar-23","2020-Mar-30")),
                      break_dates=NULL,
                      base_params,
                      data,
                      opt_pars=list(params=c(log_E0=4,
                                             log_beta0=-1),
                                    logit_rel_beta0=c(-1,-1),
                                    log_nb_disp=NULL),
                      fixed_pars=NULL,
                      sim_fun=run_sim_break,
                      sim_args=NULL,
                      aggregate_args=NULL,
                      condense_args=NULL,
                      priors=NULL,
                      debug=FALSE,
                      debug_plot=FALSE,
                      debug_hist=FALSE,
                      last_debug_plot=FALSE,
                      use_DEoptim=FALSE,
                      mle2_method="Nelder-Mead",
                      mle2_control=list(maxit=10000),
                      mle2_args=list(),
                      seed=NULL,
                      DE_args=list(),
                      DE_lwr=NULL,
                      DE_upr=NULL,
                      DE_cores=getOption("mc.cores",2))
{
    start_time <- proc.time()
    if (!is.null(break_dates)) {
        warning("use of break_dates as a top-level parameter is deprecated: please use time_args=list(break_dates=...)")
        time_args <- list(break_dates=break_dates)
    }
    v <- na.omit(data$value)
    if (any(abs(v-round(v))>1e-9)) {
        stop("need integer values in reported data (to match dnbinom)")
    }
    cc <- match.call()
    if (debug) {
        cat("start date: ", format(start_date),
            "; end_date: ", format(end_date), "\n")
    }
    ## FIXME: check that appropriate var names are present
    ## translate state variables names in data to expected internal names (e.g. newConfirmations -> reports)
    ## FIXME: should this be external?
    data$var <- trans_state_vars(data$var)
    ## FIXME: test this clause ...
    if (is.null(opt_pars$log_nb_disp) && !any(grepl("^obs_disp", names(base_params)))) {
        nvar <- length(var_names <- unique(data$var))
        opt_pars$log_nb_disp <- setNames(rep(0,nvar),var_names)
    }
    value <- pred <- NULL ## global var check
    ## FIXME: check order of dates?
    ## FIXME: allow for aggregation of reports etc.
    ## MLE
    mle_args <- nlist(start_date
                    , end_date
                    , time_args
                    , base_params
                    , opt_pars
                    , sim_args
                    , aggregate_args
                    , debug
                    , debug_hist
                    , debug_plot
                      ## FIXME: refactor debugging args!
                    , data
                    , priors
                    , sim_fun)
    opt_inputs <- unlist(opt_pars)


    debug_env$opt_inputs <- opt_inputs
    ## initialize debug history
    if (debug_hist) {
        with(debug_env,  {
            debug_hist_chunksize <- 1e4
            base_mat <- history_mat <- matrix(NA,
                                              nrow=debug_hist_chunksize,
                                              ncol=length(opt_inputs)+1,
                                              dimnames=list(NULL,c(names(opt_inputs),"NLL")))
            debug_ctr <- 1
        }) ## with()
    }        
    
    de_cal1 <- de_time <- NULL
    if (use_DEoptim) {
        DE_lims <- get_DE_lims(opt_pars)
        if (is.null(DE_lwr)) {
            DE_lwr <- DE_lims[["lwr"]]
        }
        if (is.null(DE_upr)) {
            DE_upr <- DE_lims[["upr"]]
        }
        de_ctrl_args <- list(packages=list("McMasterPandemic","bbmle"))
        if (!is.null(DE_args$control)) {
            for (n in names(DE_args$control)) {
                de_ctrl_args[[n]] <- DE_args$control[[n]]
            }
            DE_args$control <- NULL
        }
        ## test
        ## do.call(mle_fun, c(list(opt_inputs), mle_args))
        de_ctrl <- do.call(DEoptim::DEoptim.control,de_ctrl_args)
        de_arglist <- c(list(fn=mle_fun
                           , lower=DE_lwr
                           , upper=DE_upr
                           , control=de_ctrl)
                      , DE_args
                      , mle_args)
        if (DE_cores>1) {
            cl <- parallel::makeCluster(DE_cores)
            de_arglist$control$cluster <- cl
        }
        old_opts <- options(MP_badsum_action="ignore")
        de_time <- system.time(de_cal1 <- do.call(DEoptim::DEoptim,de_arglist))
        options(old_opts)
        if (DE_cores>1) parallel::stopCluster(cl)
        opt_inputs <- de_cal1$optim$bestmem   ## set input for MLE to best
        M <- de_cal1$member$pop   ## find the population of parameter sets at the last step
        ## FIXME: can we avoid re-running likelihoods to threshold?
        nll_vals <- NULL
        ## recalculate neg log-likelihoods
        nll_vals <- apply(M,1,
                          function(x) do.call(mle_fun,c(list(x),mle_args)))
        ## weighted covariance based on *likelihoods*
        likvec <- exp(-(nll_vals-min(nll_vals)))
        vM <- stats::cov.wt(M,wt=pmax(likvec,min(likvec[likvec>0])))$cov
        dimnames(vM) <- list(names(opt_inputs), names(opt_inputs))
        ## attach to de object
        de_cal1$member$Sigma <- vM
        de_cal1$member$nll_vals <- nll_vals
    } ## end of DEoptim loop
    ## now 'polish' with mle2
    ## n.b.: this has to be hacked dynamically ...
    ## FIXME: same as mle_args?
    parnames(mle_fun) <- names(opt_inputs)
    mle_data <- mle_args
    ## FIXME: is this even necessary?
    mle_data$fixed_pars <- mle_data$fixed_pars
    mle2_args <- c(mle2_args,list(namedrop_args=FALSE))
    opt_args <- c(list(minuslogl=mle_fun
                        , start=opt_inputs
                        , data=mle_data
                        , vecpar=TRUE
                        , method=mle2_method
                        , control=mle2_control
                       ),
                  mle2_args)
    mle1 <- do.call(mle_fun,c(list(opt_inputs), mle_args))
    mle2_time <- system.time(opt <- do.call(bbmle::mle2,opt_args))
    if (last_debug_plot) {
        grDevices::pdf(".debug_plot.pdf")
        mle_args$debug_plot <- TRUE
        do.call(mle_fun,c(list(coef(opt)), mle_args))
        grDevices::dev.off()
    }
    res <- list(mle2=opt
              , forecast_args = mle_data[setdiff(names(mle_data),
                            ## leave these out, they're needed for mle but not forecast
                            c("debug","debug_plot","fixed_pars","data","priors"))]
              , call = cc)
    end_time <- proc.time()
    if (debug_hist) {
        with(debug_env,
             history_mat <- history_mat[rowSums(is.na(history_mat))<ncol(history_mat),]
        )
        attr(res,"debug_hist") <- debug_env$history_mat
    }
    attr(res,"de") <- de_cal1
    attr(res,"de_time") <- de_time
    attr(res,"mle2_time") <- mle2_time
    attr(res,"total_time") <- end_time-start_time
    class(res) <- "fit_pansim"
    return(res)
}

##' find confidence envelopes by simulation
##' @inheritParams forecast_sim
##' @param fit output from \code{calibrate}
##' @param nsim number of simulations
##' @param seed random-number seed
##' @param forecast_args arguments to pass to \code{forecast_sim}
##' @param imp_wts use importance weighting, i.e. weight ensemble based on log-likelihood?
##' @param qvec vector of quantiles: NULL to return an array (nt x nvars x nsims) instead of a tibble with date/var+ quantiles
##' @param raw_ensembles (logical) return ensembles (FIXME: should implement return_type=c("array","quantiles","raw") (but not really raw, "dataframe"??))
##' @param qnames quantile names
##' @param fix_pars_re a regular expression specifying the names of parameters that should be treated as fixed when constructing the parameter ensemble
##' @param .progress progress bar?
##' @param Sigma covariance matrix to pass to \code{pop_pred_samp}
##' @param scale_Sigma multiplier for covariance matrix
##' @export
## FIXME: way to add args to forecast_args list, e.g. stochastic components?
forecast_ensemble <- function(fit,
                              nsim=200,
                              forecast_args=fit$forecast_args,
                              qvec=c(0.05,0.5,0.95),
                              qnames=c("lwr","value","upr"),
                              seed=NULL,
                              imp_wts=FALSE,
                              Sigma=bbmle::vcov(fit$mle2),
                              scale_Sigma=1,
                              calc_Rt=FALSE,
                              fix_pars_re="nb_disp",
                              raw_ensembles = FALSE,
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
                c(nlist(p,return_val,calc_Rt),forecast_args))
    }

    ## baseline fit
    r <- ff(coef(fit$mle2), return_val="aggsim")

    ## FIXME: count number of distribution params?


    if (is.null(fix_pars_re)) {
        fix_pars <- NULL
    } else {
        fix_pars <- grep(fix_pars_re, names(coef(fit$mle2)), value=TRUE)
    }

    f_args <- forecast_args
    f_args <- f_args[!names(f_args) %in% c("stoch", "stoch_start", "fixed_pars")]
    pps_args <- c(nlist(object=fit$mle2
                     , n=nsim
                     , PDify =TRUE
                     , Sigma
                     , return_wts=imp_wts
                     , data=fit$mle2@data$data
                     , fix_params=fix_pars
                     , scale_Sigma
                       )
                , f_args)


    ## use internal pop_pred_samp not bbmle version
    e_pars <- do.call(pop_pred_samp, pps_args)

    wts <- rep(1, nsim)
    if (imp_wts) {
        wts <- e_pars[,"wts"]
        if ((es <- attr(e_pars,"eff_samp"))<10) warning("low effective sample size of importance weights",
                                                        sprintf(" (sample=%d, eff sample=%1.1f)",
                                                                nrow(e_pars),es))
    }

    ## e_pars <- e_pars[,setdiff(colnames(e_pars),fix_pars)]

    ## run for all param vals in ensemble
    ## tried with purrr::pmap but too much of a headache
    t1 <- system.time(e_res <- plyr::alply(as.matrix(e_pars)
                                         , .margins=1
                                         , .fun=ff
                                         , .progress=.progress  ))

    if (is.null(qvec)) {
        ## return as array
        nv <- length(unique(r$var))
        nt <- length(unique(r$date))
        aa <- array(unlist(e_res),dim=c(nv,nt,nsim))
        dimnames(aa) <- list(var=unique(r$var),
                             date=format(unique(r$date)),
                             sim=seq(nsim))
        if (imp_wts) attr(aa,"imp_wts") <- wts
        return(aa)
    } else {
        e_res2 <- e_res %>% dplyr::bind_cols()
        if(raw_ensembles){
            e0 <- (dplyr::select(r,date,var)
                   %>% dplyr::as_tibble())
            e_res2 <- dplyr::bind_cols(e0, e_res2)
            return(e_res2)
        }

        ## get quantiles per observation
        ## safe version of wtd quantile
        wq <- function(x,weights,probs) {
            bad <- !is.finite(x)
            if (all(bad)) {
                return(rep(NA,length.out=length(probs)))
            }
            return(Hmisc::wtd.quantile(x[!bad], weights[!bad], probs))
        }
        if (imp_wts) {
            e_res3 <- apply(e_res2,1,FUN=wq,weights=wts,probs=qvec) ## na.rm is handled internally (automatically)
        } else {
            e_res3 <- apply(e_res2,1,stats::quantile,probs=qvec, na.rm=TRUE)
        }        
        e_res4 <- (e_res3 
            %>% t()
            %>% as.data.frame() ## create names (without messages); avoid need for .name_repair
            %>% dplyr::as_tibble()
            %>% setNames(qnames)
        )
    
        ## date/var values
        e0 <- (dplyr::select(r,date,var)
            %>% dplyr::as_tibble()
        )
        ## combine quantiles with the original date/var columns
        e_res3 <- dplyr::bind_cols(e0, e_res4)
        return(e_res3)
    }
}

## make a factor/dummy variables for which 'mobility period' we're in,
## i.e. before first, between first and second, etc...
## count the number of breaks in the vector that are *earlier* than the current date in X_dat
## i.e. 0, 1, 2, 3 ... make this a factor
## add it to X_dat
## d1 <- seq.Date(as.Date("2020-01-05"),as.Date("2020-01-12"),by="1 day")
## brks <- c(as.Date("2020-01-09"),as.Date("2020-01-11"))
date_count <- function(d1,brks) {
    res <- numeric(length(d1))
    for (i in seq_along(brks)) {
        res <- res + as.numeric(d1>brks[i])
    }
    return(res)
}

date_logist <- function(date_vec, date_prev, date_next=NA,
                        scale) {
    t1 <- as.numeric(date_vec-date_prev)
    r <- rep(1,length(date_vec))
    if (!is.na(date_prev)) r <- plogis(t1/scale)
    if (!is.na(date_next)) {
        t2 <- as.numeric(date_vec-date_next)
        r <- r - plogis(t2/scale)
    }
    return(r)
}

##' Combined calibration of model to multiple types of data
##'
##' Top-level calibration based on mobility, splines, and
##' phenomenological heterogeneity.  This function is a wrapper for
##' \code{\link{calibrate}}, which first sets up the model matrix for
##' time-dependent transmission rates \eqn{\beta(t)}.  For example,
##' for splines, it sets up the spline basis and the coefficients of
##' each basis function.  For mobility, it sets up the slope and
##' intercept parameters (since mobility is treated as a log-linear
##' function between break points).
##' 
##' @param params parameters
##' @param maxit maximum iterations for Nelder-Mead/optimization step
##' @param skip.hessian skip Hessian calculation?
##' @param ... extra args
##' @param mob_data mobility data
##' @param mob_breaks vector of breakpoints for piecewise mobility
##'     model
##' @param mob_breaks_int (logical) specifies whether the intercept of
##'     the mobility/transmission relationship changes at each
##'     mobility breakpoint. The default (\code{FALSE}) specifies that
##'     mobility in 'mobility period' \code{i} is
##'     \code{beta0*(rel_mobility)^p_i}; \code{TRUE} would specify
##'     \code{beta0*(a_i)*(rel_mobility)^p_i}, with \code{a_1==1}.
##' @param mob_logist_scale specifies the scale of the smooth
##'     (logistic) transition between mobility periods, in days: if it
##'     is \code{NA} (default), the model uses sharp breaks for a
##'     piecewise-constant model.  Otherwise, the parameter used makes
##'     a logistic transition between the two periods with the
##'     specified scale.
##' @param spline_days days between spline knots
##' @param spline_df overall spline degrees of freedom
##' @param spline_setback days before end of time series to set
##'     boundary knots for spline (this implies \emph{linear}
##'     extrapolation after knots if \code{spline_type="ns"} is
##'     specified, which is probably wise)
##' @param knot_quantile_var variable to use cum dist for knotspacing
##' @param spline_pen penalization for spline
##' @param spline_type spline type ("ns" for natural spline or "bs"
##'     for b-spline)
##' @param spline_int spline intercept (??)
##' @param spline_extrap spline extrapolation model ("linear" or
##'     "constant")
##' @param testing_data data frame with columns containing dates
##'     (\code{Date}) and testing intensity (\code{intensity}) (=
##'     tests per capita per day)
##' @param use_mobility include mobility as a covariate in the model?
##' @param use_phenomhet include phenomenological heterogeneity?
##' @param use_testing include variation in testing intensity?
##' @param use_spline include spline?
##' @param vars which vars to use? (default is all in data)
##' @param return_val "fit" (return calibrated value); "X"
##'     (short-circuit/return model matrix?); "formula" (return
##'     log-linear formula for time-varying beta)
##' @param start_date start date
##' @importFrom stats quantile reformulate model.matrix
##' @importFrom dplyr distinct select
##' @importFrom tidyr drop_na
##' @importFrom stats plogis
##' @importFrom splines bs
##' @examples
##' if (require(dplyr)) {
##'   dd <- ont_all %>% trans_state_vars() %>%
##'        filter(var %in% c("H","report"))
##'  params <- read_params("ICU1.csv")
##'  ## quick and dirty example (maximize speed)
##' \dontrun{
##'  calibrate_comb(data=dd, params=params,
##'                use_spline=TRUE,
##'               maxit=10, skip.hessian=TRUE, use_DEoptim =FALSE)
##' X <-  calibrate_comb(data=dd, params=params,
##'               use_spline=TRUE,
##'               spline_type="ns",
##'               spline_setback=1,
##'               spline_extrap="constant",
##'               return_val="X")
##' matplot(X, ylab="")
##' form <-  calibrate_comb(data=dd, params=params,
##'               use_spline=TRUE,
##'               spline_type="ns",
##'               spline_setback=1,
##'               spline_extrap="constant",
##'               return_val="formula")
##' print(form)
##' }
##' }
##' @inheritParams calibrate
##' @importFrom splines ns bs
##' @export
calibrate_comb <- function(data,
                     params,
                     opt_pars=NULL,
                     mob_data=NULL,
                     mob_breaks=NULL,
                     mob_breaks_int=FALSE,
                     mob_logist_scale=NA,
                     spline_days=14,
                     spline_setback=0,
                     spline_extrap=c("linear","constant"),
                     spline_df=NA,
                     knot_quantile_var=NA,
                     spline_pen=0,
                     spline_type="bs",
                     spline_int=FALSE,
                     testing_data=NULL,
                     maxit=10000,
                     skip.hessian=FALSE,
                     use_DEoptim=TRUE,
                     DE_cores=1,
                     use_mobility=FALSE,
                     use_phenomhet=FALSE,
                     use_spline=FALSE,
                     use_testing=FALSE,
                     vars=NULL,
                     debug_plot=interactive(),
                     debug=FALSE,
                     debug_hist=FALSE,
                     return_val=c("fit","X","formula","args","time_args"),
                     start_date=NULL,
                     priors=list(),
                     ...) {
    spline_extrap <- match.arg(spline_extrap)
    return_val <- match.arg(return_val)
    t_vec <- value <- NULL ## global var check
    ## choose variables
    if (!is.null(vars)) {
        indiv_vars <- trimws(unlist(strsplit(vars,"/")))
        dvars <- unique(data$var)
        if (any(!indiv_vars %in% dvars)) {
            stop("vars requested that are not in data")
        }
        data <- filter(data, var %in% indiv_vars)
    }
    if (is.null(opt_pars)) {
        opt_pars <- get_opt_pars(params,vars=unique(data$var))
    }
    if (use_phenomhet) opt_pars$params <- c(opt_pars$params,log_zeta=1)
    ## need tmp_dat ouside use_spline for non-spline, non-mob cases
    X_dat <- (data
        %>% select(date)
        %>% distinct()
        %>% mutate(t_vec=as.numeric(date-min(date)))
    )
    loglin_terms <- "-1"
    if (use_mobility) {
        if (is.null(mob_breaks)) {
            ## single activity term for the whole time period
            loglin_terms <- c(loglin_terms, "log(rel_activity)")
        } else {
            mob_breaks <- as.Date(mob_breaks)
            if (is.na(mob_logist_scale)) {
                ## everything is piecewise
                mob_breaks <- as.Date(mob_breaks)
                mob_count <- date_count(X_dat$date,mob_breaks)
                X_dat$mob_breaks <- factor(mob_count)
                ## set up changes in mobility as piecewise constant
                loglin_terms <- c(loglin_terms, "log(rel_activity):mob_breaks")
                if (mob_breaks_int) {
                    ## need to work around existing intercept term ...
                    ## create dummy variables, skip the intercept
                    tmpX <- model.matrix(~mob_breaks,data=X_dat)[,-1,drop=FALSE]
                    colnames(tmpX) <- paste0("mob_breaks",seq(ncol(tmpX)))
                    ## don't use cbind(), will coerce! 
                    X_dat <- data.frame(X_dat,tmpX)
                    loglin_terms <- c(loglin_terms,colnames(tmpX))
                }
            } else {
                ## set up changes in mobility (and intercept) as logistic transitions
                tmpX <- matrix(ncol=0,nrow=nrow(X_dat))
                m <- c(NA,mob_breaks)
                for (i in seq_along(m)) {
                    tmpX <- cbind(tmpX,date_logist(X_dat$date,
                                                   m[i],
                                                   m[i+1],
                                                   mob_logist_scale))
                }
                colnames(tmpX) <- paste0("mob_logist",seq_along(m)-1)
                loglin_terms <- c(loglin_terms,paste0("log(rel_activity):",colnames(tmpX)))
                X_dat <- data.frame(X_dat,tmpX)
                if (mob_breaks_int) {
                    ## leave out intercept
                    loglin_terms <- c(loglin_terms,colnames(tmpX)[-1])
                }
            } ## use logistic transitions
        } ## >1 period
    } ## use mobility
    if (use_spline) {
        if (is.na(spline_df)) {
            spline_df <- round(length(X_dat$t_vec)/spline_days)
        }
        knot_args <- "t_vec"
        if (is.na(knot_quantile_var)) {
            knot_args <- c(knot_args, "df=spline_df")
        } else {
            ## df specified: pick (df-degree) knots
            knot_dat <- (filter(data,var==knot_quantile_var)
                %>% na.omit()
                %>% mutate(cum=cumsum(value),t_vec=as.numeric(date-min(date)))
            )
            q_vec <- quantile(knot_dat$cum,seq(1,spline_df-3)/(spline_df-2))
            ## find closest dates to times with these cum sums
            knot_t_vec <- with(knot_dat,approx(cum,t_vec,xout=q_vec,ties=mean))$y
            ## cat("knots:",knot_t_vec,sep="\n","\n")
            knot_args <- c(knot_args, "knots=knot_t_vec")
        }
        if (spline_setback>0) {

            ## this makes the spline extrapolate linearly for the last (spline_setback) days of the time series,
            ## in a reasonably elegant way, if we use spline_type="ns"
            ##
            ## if we want to work harder and make the spline component *constant* up to (spline_setback) days
            ## before the end of the data we will have to hack more.
            ##
            ## strategy: take the generated model matrix and, for the columns corresponding to the spline basis,
            ## replace the last (spline_setback-1) rows with row (nt-spline_setback)
            ## 
            ## This has a slightly bad interaction with the specification of the number/placement of knots
            ## (i.e. we should probably allocate and place knots on the basis of the shortened time series as well),
            ##  but too horrible to think about for right now ... so we will end up with slightly fewer knots
            knot_args <- c(knot_args, "Boundary.knots=c(min(t_vec),max(t_vec)-spline_setback)")
        }
        spline_term <- sprintf("%s(%s)",spline_type,paste(knot_args,collapse=","))
        if(spline_int){
            spline_term <- c("spline_int",spline_term)
            X_dat$spline_int <- 1
        }
        loglin_terms <- c(loglin_terms, spline_term)
    }
    form <- reformulate(loglin_terms)
    if (use_mobility) {
        X_dat <- (full_join(X_dat,mob_data,by="date")
            %>% arrange(date)  ## fill_edge_values assumes ordered!
            %>% drop_na(t_vec)  ## omit values *before* data start ...
            ## FIXME: might want to keep these, but then we have to worry about how to set
            ##  spline values constant before the data start (as we don't want to waste
            ##  spline knots out there ...)
            %>% mutate_at("rel_activity",fill_edge_values)
        )
    }
    if (return_val=="formula") {
        return(form)
    }
    if (use_testing) {
        params <- update(params,testing_intensity=testing_data$intensity[1])
        ## set up in 'timevars' format (passed through run_sim_loglin to run_sim)
        testing_data <- with(testing_data,
                             data.frame(Date=Date,
                                        Symbol="testing_intensity",
                                        Relative_value=c(1, intensity[-1]/intensity[1])))
    }
    X <- model.matrix(form, data = X_dat)
    if (use_spline && spline_setback>0 && spline_extrap=="constant") {
        spline_cols <- grep(paste0(spline_type,"("), fixed=TRUE, colnames(X))
        last_t <- max(X_dat$t_vec)-spline_setback
        X[X_dat$t_vec>last_t,spline_cols] <- matrix(rep(X[X_dat$t_vec==last_t,], spline_setback),
                                                     ncol=length(spline_cols),
                                                     byrow=TRUE)
    }
    if (return_val=="X") return(X)
    if(is.null(opt_pars$time_beta)){
    ## matplot(X_dat$t_vec,X,type="l",lwd=2)
    opt_pars$time_beta <- rep(0,ncol(X))  ## mob-power is incorporated (param 1)
    names(opt_pars$time_beta) <- colnames(X)
    }
    time_args <- nlist(X,X_date=X_dat$date)
    if (use_testing) {
        time_args <- c(time_args, list(testing_data = testing_data))
	      if(return_val == "time_args")return(time_args)
    }
    if (use_spline && spline_pen>0) {
        spline_pars <- grep("^[bn]s\\(",names(opt_pars$time_beta))
        spline_beg <- spline_pars[1]
        spline_end <- spline_pars[length(spline_pars)]
        priors <- c(priors,list(bquote(~sum(dnorm(time_beta[.(spline_beg):.(spline_end)],mean=0,sd=.(1/spline_pen))))))
    }
    ## do the calibration
    ## debug <- use_DEoptim
    argList <- c(nlist(data
                     , use_DEoptim
                     , DE_cores
                     , debug_plot
                     , debug
                     , debug_hist
                     , base_params=params
                     , mle2_control = list(maxit=maxit)
                     , mle2_args=list(skip.hessian=skip.hessian)
                     , opt_pars
                     , time_args
                     , sim_fun = run_sim_loglin
                     , priors)
               , list(...))

    if (return_val=="args") return(argList)
    
    res <- do.call(calibrate, argList)

    ## FIXME: figure out descriptive string?

    ## checkpoint
    ## saveRDS(res,file=sprintf("ont_fac_%s_fit.rds",flags))
    return(res)  ## return *after* saving!!!!
}

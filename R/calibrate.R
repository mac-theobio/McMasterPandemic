
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
##' @inheritParams run_sim
##' @param p vector of parameters - on the link (log/logit) scale as appropriate; these are the parameters that will be adjusted during calibration
##' @param condense_args arguments to pass to \code{\link{condense}} (via \code{\link{run_sim}})
##' @param stoch stochastic settings (see \code{\link{run_sim}})
##' @param return_val specify values to return (aggregated simulation, or just the values?)
##' @export
forecast_sim <- function(p, opt_pars, base_params, start_date, end_date, break_dates,
                         fixed_pars = NULL,
                         stoch = NULL,
                         stoch_start = NULL,
                         sim_args=NULL, aggregate_args=NULL, condense_args=NULL,
                         ## FIXME: return_val is redundant with sim_fun
                         return_val=c("aggsim","vals_only"),
                         debug = FALSE)
{
    return_val <- match.arg(return_val)
    if (!is.null(stoch)) {
        if (is.null(sim_args)) {
            sim_args <- nlist(stoch, stoch_start)
        } else {
            sim_args <- c(sim_args, nlist(stoch=stoch, stoch_start))
        }
    }
    ## restructure and inverse-link parameters
    ## if (debug) cat("forecast 0",opt_pars[["log_beta0"]],"\n")
    pp <- invlink_trans(restore(p, opt_pars, fixed_pars))
    ## substitute into parameters
    params <- update(base_params, params=pp$params, .list=TRUE)
    ## if (debug) cat("forecast ",params[["beta0"]],"\n")
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

##' negative log-likelihood function
##' @param p parameter vector (in unlisted form)
##' @param ... unused (but useful in case junk needs to be discarded)
##' @inheritParams calibrate
##' @export
mle_fun <- function(p, data, debug=FALSE, debug_plot=FALSE,
                    opt_pars, base_params, start_date, end_date, break_dates=NULL,
                    sim_args=NULL, aggregate_args=NULL,
                    priors=NULL, ...) {
    ## ... is to drop any extra crap that gets in there
    if (debug) cat("mle_fun: ",p,"\n")
    var <- pred <- value <- NULL    ## defeat global-variable checkers
    ## pass everything to forecaster
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
            with(subset(r2,value>0),
                 plot(date,value,col=as.numeric(factor(var)),log="y"))
            r2s <- split(r2,r2$var)
            invisible(Map(function(x,c) lines(x$date,x$pred, col=c), r2s, seq_along(r2s)))
        }
    }
    ## need to reconstruct parameter set to get NB parameter(s)
    ## FIXME: fixed params can now be handled through mle2?
    pp <- invlink_trans(restore(p, opt_pars))
    if (length(names(pp$nb_disp))>0) {
        ## FIXME: invlink_trans doesn't handle this case correctly!
        ##  for now, manually exp() ...
        pp$nb_disp <- exp(pp$nb_disp)
    }
    if (length(pp$nb_disp)==1) {
        dvals <- with(r2,dnbinom(value,mu=pred,size=pp$nb_disp,log=TRUE))        
    } else {
        r2 <- merge(r2,data.frame(var=names(pp$nb_disp),nb_disp=pp$nb_disp),
                    by="var")
        dvals <- with(r2,dnbinom(value,mu=pred,size=exp(nb_disp),log=TRUE))
    }
    ## clamp NaN/NA values to worst obs
    dvals[is.na(dvals)] <- min(dvals,na.rm=TRUE)
    ret <- -sum(dvals)
    if (!is.null(priors)) {
        for (pr in priors) {
            pr <- pr[[2]] ## drop tilde
            pr$log <- TRUE
            ret <- ret - eval(pr,list2env(pp))
        }
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
##' @param base_params baseline parameters (an object (vector?) of type \code{params_pansim} containing all of the parameters needed for a simulation; some may be overwritten during the calibration process)
##' @param data a data set to compare to, containing date/var/value (current version assumes that only a single state var is included)
##' @param opt_pars starting parameters (and structure).  Parameters that are part of the \code{params_pansim} parameter vector can be specified within the \code{params} element (with prefixes if they are transformed); other parameters can include distributional parameters or time-varying parameters
##' @param fixed_pars parameters to fix
##' @param sim_args additional arguments to pass to \code{\link{run_sim}}
##' @param aggregate_args arguments passed to \code{\link{aggregate.pansim}}
##' @param mle2_control control args for mle2
##' @param mle2_method method arg for mle2
##' @param mle2_args additional arguments for mle2
##' @param debug print debugging messages?
##' @param debug_plot plot debugging curves? (doesn't work with parallel DEoptim)
##' @param priors a list of tilde-delimited expressions giving prior distributions expressed in terms of the elements of \code{opt_pars}, e.g. \code{list(~dlnorm(rel_beta0[1],meanlog=-1,sd=0.5))}
##' @param seed random-number seed (for DE)
##' @param use_DEoptim use differential evolution as first stage?
##' @param DE_args arguments for \code{\link{DEoptim}}
##' @param DE_lwr lower bounds for DE optimization
##' @param DE_upr upper bounds, ditto
##' @param DE_cores number of parallel workers for DE
##' @importFrom graphics lines
##' @importFrom bbmle parnames<- mle2
## DON'T import stats::coef !
##' @examples
##' library(dplyr)
##' params <- fix_pars(read_params("ICU1.csv"))
##' opt_pars <- list(params=c(log_E0=4, log_beta0=-1,
##'           log_mu=log(params[["mu"]]), logit_phi1=qlogis(params[["phi1"]])),
##'                                   logit_rel_beta0=c(-1,-1),
##'                                   log_nb_disp=NULL)
##' dd <- (ont_all %>% trans_state_vars() %>% filter(var %in% c("report", "death", "H")))
##' \dontrun{
##'    cal1 <- calibrate(data=dd, base_params=params, opt_pars=opt_pars, debug_plot=TRUE)
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
                      break_dates=c("2020-Mar-23","2020-Mar-30"),
                      base_params,
                      data,
                      opt_pars=list(params=c(log_E0=4,
                                             log_beta0=-1),
                                    logit_rel_beta0=c(-1,-1),
                                    log_nb_disp=NULL),
                      fixed_pars=NULL,
                      sim_args=NULL,
                      aggregate_args=NULL,
                      priors=NULL,
                      debug=FALSE,
                      debug_plot=FALSE,
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
    cc <- match.call()
    if (debug) {
        cat("start date: ", format(start_date),
            "; end_date: ", format(end_date), "\n")
    }
    ## FIXME: check that appropriate var names are present
    ## translate state variables names in data to expected internal names (e.g. newConfirmations -> reports)
    ## FIXME: should this be external?
    data$var <- trans_state_vars(data$var)
    if (is.null(opt_pars$log_nb_disp)) {
        nvar <- length(var_names <- unique(data$var))
        opt_pars$log_nb_disp <- setNames(rep(0,nvar),var_names)
    }
    value <- pred <- NULL ## global var check
    ## FIXME: check order of dates?
    ## FIXME: allow for aggregation of reports etc.
    ## MLE
    mle_args <- nlist(start_date, end_date,
                      break_dates,
                      base_params,
                      opt_pars,
                      sim_args,
                      aggregate_args,
                      debug,
                      debug_plot,
                      data,
                      priors)
    opt_inputs <- unlist(opt_pars)
    de_cal1 <- NULL
    if (use_DEoptim) {
        if (is.null(DE_lwr)) {
            DE_lwr <- opt_inputs
            DE_lwr[] <- -1
            DE_lwr[["params.log_E0"]] <- 1
        }
        if (is.null(DE_upr)) {
            DE_upr <- opt_inputs
            DE_upr[] <- 1
            DE_upr[grepl("rel_beta0",names(DE_upr))] <- 4
            DE_upr[grepl("nb_disp|E0",names(DE_upr))] <- 5
        }
        cl <- parallel::makeCluster(DE_cores)
        de_arglist <- c(list(fn=mle_fun, lower=DE_lwr, upper=DE_upr,
                             control=DEoptim::DEoptim.control(storepopfrom=1,
                                                              cluster=cl,
                                                packages=list("McMasterPandemic","bbmle"))),
                        DE_args,mle_args)
        de_time <- system.time(de_cal1 <- do.call(DEoptim::DEoptim,de_arglist))
        parallel::stopCluster(cl)
        opt_inputs <- de_cal1$optim$bestmem   ## set input for MLE to best
        M <- de_cal1$member$pop
        ## FIXME: can we avoid re-running likelihoods to threshold?
        nll_vals <- apply(M,1,
                          function(x) do.call(mle_fun,c(list(x),mle_args)))
        M <- M[log10(nll_vals-min(nll_vals))<1.5,]
        ## attach to de object
        de_cal1$member$Sigma <- var(M)
    }
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
          data,
          priors)
    ## unnecessary
    ## if (utils::packageVersion("bbmle")<"1.0.23.2") stop("please remotes::install('bbolker/bbmle') to get newest version")
    opt_args <- c(list(minuslogl=mle_fun
                        , start=opt_inputs
                        , data=mle_data
                        , vecpar=TRUE
                        , method=mle2_method
                        , control=mle2_control
                        ## , namedrop_hack = FALSE  ## unnecessary
                       ),
                  mle2_args)
    opt <- do.call(bbmle::mle2,opt_args)
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
    attr(res,"de") <- de_cal1
    class(res) <- "fit_pansim"
    return(res)
}

##' find confidence envelopes by simulation
##' @param fit output from \code{calibrate}
##' @param nsim number of simulations
##' @param seed random-number seed
##' @param forecast_args arguments to pass to \code{forecast_sim}
##' @param imp_wts use importance weighting, i.e. weight ensemble based on log-likelihood?
##' @param qvec vector of quantiles
##' @param qnames quantile names
##' @param fix_pars_re a regular expression specifying the names of parameters that should be treated as fixed when constructing the parameter ensemble
##' @param .progress progress bar?
##' @param Sigma covariance matrix to pass to \code{\link{pop_pred_samp}}
##' @export
## FIXME: way to add args to forecast_args list, e.g. stochastic components?
## FIXME: use bbmle::pop_pred_samp?
forecast_ensemble <- function(fit,
                              nsim=200,
                              forecast_args=fit$forecast_args,
                              qvec=c(0.05,0.5,0.95),
                              qnames=c("lwr","value","upr"),
                              seed=NULL,
                              imp_wts=FALSE,
                              Sigma=bbmle::vcov(fit$mle2),
                              fix_pars_re="nb_disp",
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

    ## FIXME: count number of distribution params?


    if (is.null(fix_pars_re)) {
        fix_pars <- NULL
    } else {
        fix_pars <- grep(fix_pars_re, names(coef(fit$mle2)), value=TRUE)
    }

    pps_args <- c(list(fit$mle2
                     , n=nsim
                     , PDify =TRUE
                     , Sigma = Sigma
                     , return_wts=imp_wts
                     , data=fit$mle2@data$data
                     , fix_params=fix_pars)
                , forecast_args[c("opt_pars","base_params","start_date","end_date")])


    e_pars <- do.call(bbmle::pop_pred_samp, pps_args)

    wts <- rep(1, nsim)
    if (imp_wts) {
        wts <- e_pars[,"wts"]
        if (attr(e_pars,"eff_samp")<10) warning("low effective sample size of importance weights")
    }

    np <- length(unlist(forecast_args$opt_pars[c("params","logit_rel_beta0")]))
    e_pars <- e_pars[,seq(np)]

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

        ## get quantiles per observation
        if (imp_wts) {
            e_res3 <- apply(e_res2,1,Hmisc::wtd.quantile,weights=wts,probs=qvec,na.rm=TRUE)
        } else {
            e_res3 <- apply(e_res2,1,stats::quantile,probs=qvec,na.rm=TRUE)
        }        
        e_res4 <- (e_res3 
            %>% t()
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

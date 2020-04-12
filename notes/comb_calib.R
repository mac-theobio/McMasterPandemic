## this is a generalization/adaptation of 'get_break'; should
## think more generally about the interface, but for now just using
## this for the special case of RSA forecasting (fit E0, beta0,
##  breakpoints).  Use Ontario data as example.


forecast_sim <- function(p, opt_pars, base_params, start_date, end_date, break_dates,
                         sim_args=NULL, aggregate_args=NULL,
                         ## FIXME: return_val is redundant with sim_fun
                         return_val=c("aggsim","vals_only"))
{
    return_val <- match.arg(return_val)
    ## restructure and inverse-link parameters
    pp <- invlink_trans(relist(p, opt_pars))
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

## EXAMPLES: see below
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
        pp <- invlink_trans(relist(p, opt_pars))
        ret <- with(r2,-sum(dnbinom(value,mu=pred,size=pp$nb_disp,log=TRUE)))
        ## FIXME: add evaluation number?
        if (debug) cat(unlist(pp),ret,"\n")
        return(ret)
    }
    ## use optim to start with; maybe switch to mle2 later
    do.call(optim,
            c(list(par=unlist(opt_pars), fn=mle_fun, data=data,
                   do_plot=debug_plot),
              optim_args))
}

run_stuff <- TRUE
use_hosp <- FALSE
weekly <- TRUE

if (run_stuff) {
    library(McMasterPandemic)
    library(ggplot2); theme_set(theme_bw())
    library(dplyr)
    source("notes/ontario_clean.R") ## n.b. need to fix expectation of working directory; add an ON data set to pkg?
    dd <- dplyr::filter(ont_recent,var==if (!use_hosp) "newConfirmations" else "Hospitalization")
    if(weekly){
      dd <- (dd 
        %>% mutate(week = format(date, "%Y-%U"))
        %>% group_by(week)
        %>% summarise(date = max(date)
            , value = sum(value)
            , var = first(var)
            )
      )
      agg_list <- list(t_agg_start="07-Mar-2020",t_agg_period="7 days",t_agg_fun=sum)
    }
    print(ggplot(dd,aes(date,value)) + geom_point() + scale_y_log10())
    ## adjust parameters to sensible generation interval
    params <- fix_pars(read_params("ICU1.csv"), target=c(Gbar=6),u_interval=c(-1,1),
                       pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a")))
    params[["N"]] <- 19.5e6  ## reset pop to Ontario
    summary(params)
    g1 <- get_break_gen(data=dd, base_params=params, debug=TRUE,
                        optim_args=list(control=list(maxit=10000),hessian=TRUE),
                        var=if (!use_hosp) "report" else "H",
                        aggregate_args = agg_list)
    ## check standard deviations
    sqrt(diag(solve(g1$hessian)))
    ## for structure
    opt_pars <- list(log_E0=4, log_beta0=-1, log_rel_beta0=c(-1,-1), log_nb_disp=0)
    ## get parameters back onto original scale
    pp <- invlink_trans(relist(g1$par, opt_pars))
    print(pp)
    bd <- ldmy(c("23-Mar-2020","30-Mar-2020"))
    r <- forecast_sim(g1$par, opt_pars,
                      base_params=params,
                      start_date=min(dd$date)-15,
                      end_date="1-Aug-2020",
                      break_dates=bd)
    ## aggregate_args = agg_list)
    ## FIXME: r can't use plot.pansim method ATM
    print(ggplot(r,aes(date,value,colour=var))
          +geom_line()
          + scale_y_log10()
          + geom_point(data=dplyr::mutate_at(dd,"var",trans_state_vars))
          + geom_vline(xintercept=bd,lty=2)
          )

    ## parameter ensemble
    set.seed(101)
    e_pars <- as.data.frame(MASS::mvrnorm(200,
                            mu=g1$par,
                            Sigma=solve(g1$hessian)))
    ## tried with purrr::pmap but too much of a headache
    t1 <- system.time(e_res <- plyr::alply(as.matrix(e_pars)
                                         , .margins=1
                                         , .fun=forecast_sim
                                         , base_params=params
                                         , start_date=min(dd$date)-15,
                                         , end_date="1-Aug-2020"
                                         , return_val="vals_only"
                                         , break_dates=bd
                                         , opt_pars = opt_pars))

    e_res2 <- (e_res %>% bind_cols()
        %>% apply(1,quantile,c(0.1,0.5,0.9),na.rm=TRUE)
        %>% t()
        %>% as_tibble()
        %>% setNames(c("lwr","value","upr"))
    )
    e0 <- dplyr::select(r,date,var) %>% as_tibble()
    ## combine quantiles with the original date/var columns
    e_res3 <- bind_cols(e0, e_res2)
    print(ggplot(e_res3, aes(date,value,colour=var,fill=var))
          + geom_line()
          + geom_ribbon(colour=NA, alpha=0.2, aes(ymin=lwr, ymax=upr))
          + geom_point(data=dplyr::mutate_at(dd,"var",trans_state_vars))
          + geom_vline(xintercept=bd,lty=2)
          + scale_y_log10(limits=c(1,NA), oob=scales::squish)
          )
}

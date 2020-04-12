## this is a generalization/adaptation of 'get_break'; should
## think more generally about the interface, but for now just using
## this for the special case of RSA forecasting (fit E0, beta0,
##  breakpoints).  Use Ontario data as example.

##' @examples
##' g1 <- get_break_gen(trace_mle=TRUE,optim_args=list(control=list(maxit=10000)))
##' (ip1 <- invlink_trans(g1$par))
##' params <- read_params("ICU1.csv")
##' summary(update(params,beta0=ip1[["beta0"]]))
##' summary(update(params,beta0=ip1[["beta0"]]*ip1[["rel_beta01"]]))
get_break_gen <- function(date0=ldmy("1-Mar-2020"),
                          end_date=ldmy("8-Apr-2020"),
                          break_dates=c("20-Mar-2020","23-Mar-2020"),
                          base_params=read_params("ICU1.csv"),
                          data=dplyr::filter(ont_recent,var=="newConfirmations"),
                          var="report",
                          opt_pars=list(log_E0=4,
                                        log_beta0=-1,
                                        log_rel_beta0=c(-1,-1),
                                        log_nb_disp=0),
                          sim_args=NULL,
                          optim_args=NULL,
                          trace_mle=FALSE)
                          
{
    ## FIXME: check that appropriate var names are present
    data$var <- trans_state_vars(data$var)
    value <- pred <- NULL ## global var check
    ## FIXME: check order of dates?
    ## brute force
    ## we are changing beta0 at a pre-specified set of break dates
    ## (by an unknown amount)
    mle_fun <- function(p,data) {
        ## restructure and inverse-link parameters
        pp <- invlink_trans(relist(p, opt_pars))
        params <- update(base_params, E0=pp[["E0"]], beta0=pp[["beta0"]])
        s0 <- make_state(params=params)
        r <- do.call(run_sim_break,
                     c(nlist(params,
                             start_date=date0,
                             end_date,
                             break_dates,
                             rel_beta0=pp$rel_beta0),
                       sim_args))
        ## run simulation, aggregate
        r <- (aggregate(r, pivot=TRUE)
            %>% dplyr::rename(pred="value")
        )
        ## match up sim results with specified data
        names(data) <- tolower(names(data)) ## ugh
        data2 <- dplyr::filter(data,var %in% unique(r$var))
        r2 <- (dplyr::left_join(data2,r,by=c("date","var"))
            %>% tidyr::drop_na(value,pred))
        ## FIXME: why do we have an NA in pred??
        ## compute negative log-likelihood
        ## FIXME assuming a single nb_disp for now
        ret <- with(r2,-sum(dnbinom(value,mu=pred,size=pp$nb_disp,log=TRUE)))
        if (trace_mle) cat(unlist(pp),ret,"\n")
        return(ret)
    }
    ## use optim to start with; maybe switch to mle2 later
    do.call(optim,
            c(list(par=unlist(opt_pars), fn=mle_fun, data=data),
              optim_args))
}

load("ontario_clean.RData")

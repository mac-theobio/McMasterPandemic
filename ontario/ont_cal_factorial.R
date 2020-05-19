library(McMasterPandemic)
library(splines)
library(dplyr)
library(parallel)

## SHORTCUT for debugging
## mle2_control <- list(maxit = 500)
mle2_control <- list(maxit=10000)

load("ontario_clean.RData")
load("ont_cal.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)


## select the part of the mobility data that lies within the calibration data set
comb_sub2 <- (comb_sub
    %>% right_join(ont_all %>% select(date) %>% unique(),by="date")
    %>% na.omit()
    %>% mutate_at("rel_activity",pmin,1) ## cap relative mobility at 1
)

print(unique(ont_all_sub$var))
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")

X_date <- unique(ont_noICU$date)
comb_sub3 <- (full_join(unique(select(ont_noICU,date)),comb_sub2,by="date")
    ## assume missing data = constant at last value
    %>% mutate_at("rel_activity",zoo::na.locf)
    %>% mutate(t_vec=as.numeric(date-min(date)))
)

# plot(rel_activity~date,data=comb_sub3)

opt_pars <- list(
    ## these params are part of the main parameter vector: go to run_sim()
    params=c(log_E0=4      ## initial exposed
           , log_beta0=-1  ## initial baseline transmission
             ## fraction of mild (non-hosp) cases
           , log_mu=log(params[["mu"]])
             ## fraction of incidence reported
             ## logit_c_prop=qlogis(params[["c_prop"]]),
             ## fraction of hosp to acute (non-ICU)
           , logit_phi1=qlogis(params[["phi1"]])
             )
  , log_nb_disp=NULL
)


run_cali <- function(flags) {

    print(flags)

    ss <- function(flags,i) as.logical(as.numeric(substr(flags,i,i)))
    ## flags
    use_DEoptim <-ss(flags,1)
    use_mobility <- ss(flags,2)
    use_spline <- ss(flags,3)
    use_zeta <- ss(flags,4)

    ## now modify opt_pars according to flags
    if (use_zeta) opt_pars$params <- c(opt_pars$params,log_zeta=1)
    loglin_terms <- "-1"
    if (use_mobility) loglin_terms <- c(loglin_terms, "log(rel_activity)")
    if (use_spline) loglin_terms <- c(loglin_terms, "bs(t_vec,df=round(length(t_vec)/14))")
    form <- reformulate(loglin_terms)
    X <- model.matrix(form, data = comb_sub3)

    matplot(comb_sub3$t_vec,X,type="l",lwd=2)
    opt_pars$time_beta <- rep(0,ncol(X))  ## mob-power is incorporated (param 1)
    ## opt_pars$time_beta <- (0:(ncol(X)-1))/10  ## hack so we can see what's going on
    time_args <- nlist(X,X_date=comb_sub3$date)

    r0 <- run_sim_loglin(params=params, extra_pars=opt_pars["time_beta"] ## n.b. single brackets here are on purpose
                 , time_args=time_args
                 , sim_args=list(start_date=min(ont_all_sub$date)-15
                               , end_date=max(ont_all_sub$date)))

    plot(r0,break_dates=NULL,log=TRUE)
    ## do the calibration
    debug <- !use_DEoptim

    t_ont_cal_fac <- system.time(ont_cal_fac <-
                                     calibrate(data=ont_noICU
                                             , use_DEoptim=use_DEoptim
                                             , DE_cores = 1
                                             , debug_plot = interactive()
                                             , debug = debug
                                             , base_params=params
                                             , mle2_control = mle2_control
                                             , opt_pars = opt_pars
                                             , time_args=time_args
                                             , sim_fun = run_sim_loglin
                                               )
                                 ) ## system.time

    
    res <- list(mod = flags, time = t_ont_cal_fac, fit = ont_cal_fac, data=ont_noICU)
    ## checkpoint
    saveRDS(res,file=sprintf("ont_fac_%s_fit.rds",flags))
    return(res)  ## return after saving!!!!
}


## testing
if (FALSE) {
    r1 <- run_cali("1001")
    r2 <- run_cali("0001")
    library(bbmle)
    r1 <- run_cali("0001")
    plot(r1$fit, data=ont_noICU)
    r2 <- run_cali("0010")
    saveRDS(r2,file="r2_tmp.rds")
    r2H <- r2
    names(r2H$fit$forecast_args$opt_pars)
    r2H$fit$forecast_args$opt_pars$log_nb_disp <- NULL
    r2H$fit$forecast_args$opt_pars
    plot(r2$fit,data=r2$data)
    plot(r2H$fit,data=r2$data)
    p <- coef(r2$fit$mle2)
    invlink_trans(restore(p,r2$fit$forecast_args$opt_pars))
    invlink_trans(restore(p,r2H$fit$forecast_args$opt_pars))
    coef(r2$fit$mle2)
    opt_pars2 <- opt_pars
    opt_pars2$time_beta <- rep(0,ncol(X))  ## mob-power is incorporated (param 1)
    ## what is different e.g. between last call to 
    f_args <- r2$fit$forecast_args
    f_args2 <- readRDS(".mle_checkpoint.rds")
    setdiff(names(f_args), names(f_args2))
    setdiff(names(f_args2), names(f_args))
    nm <- intersect(names(f_args),names(f_args2))
    setNames(lapply(nm,function(n)all.equal(f_args[[n]], f_args2[[n]])), nm)
    ff <- (do.call(forecast_sim, f_args2)
        %>% filter(var %in% McMasterPandemic:::keep_vars)
    )
    ff_pred <- (do.call(forecast_sim,
                        c(list(p=f_args2$p)
                               ## coef(r2$fit$mle2))
                             , f_args))
        %>% filter(var %in% McMasterPandemic:::keep_vars)
    )
    f_args3 <- f_args
    f_args3$opt_pars <- f_args2$opt_pars
    f_args3$p <- f_args2$p
    ff_hybrid <- (do.call(forecast_sim, f_args3)
        %>% filter(var %in% McMasterPandemic:::keep_vars)
    )
    
    print(ggplot(bind_rows(forecast=ff,pred=ff_pred,hybrid=ff_hybrid,.id="source"),
           aes(date,value,colour=var,lty=source,shape=source))
          + geom_line()
          + geom_point()
          + scale_y_log10()
          + facet_wrap(~source)
          )
    pp <- predict(r2$fit)
    r2fix <- r2
    r2fix$fit$forecast_args$opt_pars <- f_args2$opt_pars
    plot(r2fix$fit, data=ont_noICU)
}

factorial_combos <- apply(expand.grid(replicate(4,0:1,simplify=FALSE)),
                1,paste,collapse="")
res_list <- mclapply(factorial_combos, run_cali, mc.cores=5)

## run_cali("1111")
#rdsave(res_list, factorial_combos)

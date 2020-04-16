use_hosp <- FALSE
weekly <- TRUE

library(McMasterPandemic)
library(anytime)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
source("../ontario/ontario_clean.R") ## add an ON data set to pkg?
dd <- dplyr::filter(ont_recent,
          var==if (!use_hosp) "newConfirmations" else "Hospitalization")
if (weekly) {
    dd <- (ont_all 
        %>% mutate(week = format(date, "%Y-%U"))
        %>% group_by(week,var)
        %>% summarise(date = max(date)
                    , value = sum(value)
                      )
        %>% ungroup()
        %>% mutate(diff=diff(c(0,value)))
        %>% dplyr::filter(diff>0)
    )
    ## set start values to (initial date in data set - 6 days) to make *end* of first
    ## aggregation period line up correctly
    agg_list <- list(t_agg_start=min(dd$date)-6,t_agg_period="7 days",t_agg_fun=sum)
} else {
    dd <- dd %>% dplyr::filter(date>as.Date("2020-03-15"))
    agg_list <- NULL
}
print(ggplot(dd,aes(date,value)) + geom_point() + scale_y_log10())
## adjust parameters to sensible generation interval
params <- fix_pars(read_params("ICU1.csv"), target=c(Gbar=6),
                   pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a")))
params[["N"]] <- 19.5e6  ## reset pop to Ontario
summary(params)
g1 <- calibrate(data=dd, base_params=params,
                optim_args=list(control=list(maxit=10000),hessian=TRUE),
                aggregate_args = agg_list)
fp <- list(log_nb_disp=2)
g1R <- calibrate(data=dd, base_params=params,
                 optim_args=list(control=list(maxit=10000),hessian=TRUE),
                 fixed_pars=fp,
                 aggregate_args = agg_list)

## check standard deviations
sqrt(diag(solve(g1$hessian)))



## re-run forecast with best-fit parameters
f_args <-attr(g1,"forecast_args")
r <- do.call(forecast_sim,
    c(list(p=g1$par), f_args))

## FIXME: r can't use plot.pansim method ATM
keep_vars <- c("H","ICU","D","report")
rs <- dplyr::filter(r, var %in% keep_vars)
print(ggplot(rs,aes(date,value,colour=var))
      + geom_line()
      + scale_y_log10(limits=c(1,NA),oob=scales::squish)
      + geom_point(data=dplyr::mutate_at(dd,"var",trans_state_vars))
      + geom_vline(xintercept=anydate(f_args$break_dates),lty=2)
      )

e_res3 <- forecast_ensemble(g1) %>% filter(var %in% keep_vars)

print(ggplot(e_res3, aes(date,value,colour=var,fill=var))
      + geom_line()
      + geom_ribbon(colour=NA, alpha=0.2, aes(ymin=lwr, ymax=upr))
      + geom_point(data=dplyr::mutate_at(dd,"var",trans_state_vars))
      + geom_vline(xintercept=anydate(f_args$break_dates),lty=2)
      + scale_y_log10(limits=c(1,NA), oob=scales::squish)
      )

### now try calibrating to
ont_recent_t <- mutate_at(ont_recent,"var",trans_state_vars)
dd2 <- dplyr::filter(ont_recent_t, var != "Ventilator")

## fit everything: cases, hosp, ICU, death, along with params
##    specifying fractions determining those flows
##  FIXME: allow variable-specific NB disp
opt_pars2 <- list(params=c(log_E0=4, log_beta0=-1,
                           ## fraction of mild (non-hosp) cases
                           log_mu=log(params[["mu"]]),
                           ## fraction of incidence reported
                           logit_c_prop=qlogis(params[["c_prop"]]),
                           ## fraction of hosp to acute (non-ICU)
                           logit_phi1=qlogis(params[["phi1"]]),
                           ## fraction of ICU cases dying
                           logit_phi2=qlogis(params[["phi2"]])),
                           log_rel_beta0 = c(-1,-1),
                           log_nb_disp=0)
t2 <- system.time(g2 <- calibrate(data=dd2, base_params=params,
                opt_pars=opt_pars2,
                optim_args=list(control=list(maxit=10000),hessian=TRUE))
)
                           
r2 <- do.call(forecast_sim, c(list(p=g2$par), attr(g2,"forecast_args")))
(ggplot(dplyr::filter(r2, var %in% keep_vars),
        aes(date,value,colour=var))
    + geom_line()
    + scale_y_log10()
    + geom_point(data=dd2)
    + geom_vline(xintercept=anydate(attr(g2,"forecast_args")$break_dates),lty=2)
)
i1 <- invlink_trans(restore(g2$par,opt_pars2))
params_fitted <- update(params,i1$params)
## hack
moment_params <- eval(formals(fix_pars)$pars_adj)[[2]]
out <- (describe_params(params_fitted)
    %>% mutate(type=case_when(
                   symbol %in% names(i1$params) ~ "mle-calibrated",
                   symbol %in% moment_params ~ "Gbar-calibrated",
                   TRUE ~ "assumed"),
               type=factor(type,levels=c("mle-calibrated",
                                         "Gbar-calibrated",
                                         "assumed")))
    %>% arrange(type)
)

write_csv(out, path="ont_calib.csv")

use_hosp <- FALSE
weekly <- FALSE

library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
source("notes/ontario_clean.R") ## add an ON data set to pkg?
dd <- dplyr::filter(ont_recent,
          var==if (!use_hosp) "newConfirmations" else "Hospitalization")
if (weekly) {
    dd <- (dd 
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
      + geom_vline(xintercept=ldmy(f_args$break_dates),lty=2)
      )

e_res3 <- forecast_ensemble(g1) %>% filter(var %in% keep_vars)

print(ggplot(e_res3, aes(date,value,colour=var,fill=var))
      + geom_line()
      + geom_ribbon(colour=NA, alpha=0.2, aes(ymin=lwr, ymax=upr))
      + geom_point(data=dplyr::mutate_at(dd,"var",trans_state_vars))
      + geom_vline(xintercept=ldmy(f_args$break_dates),lty=2)
      + scale_y_log10(limits=c(1,NA), oob=scales::squish)
      )

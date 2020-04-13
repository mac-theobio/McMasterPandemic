use_hosp <- FALSE
weekly <- TRUE
source("notes/calib_funs.R")

library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
source("notes/ontario_clean.R") ## n.b. need to fix expectation of working directory; add an ON data set to pkg?
dd <- dplyr::filter(ont_recent,var==if (!use_hosp) "newConfirmations" else "Hospitalization")
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
fp <- list(log_nb_disp=2)
g1R <- get_break_gen(data=dd, base_params=params, debug=TRUE,
                     optim_args=list(control=list(maxit=10000),hessian=TRUE),
                     fixed_pars=fp,
                     var=if (!use_hosp) "report" else "H",
                     aggregate_args = agg_list)

## check standard deviations
sqrt(diag(solve(g1$hessian)))
## for structure
opt_pars <- list(log_E0=4, log_beta0=-1, log_rel_beta0=c(-1,-1), log_nb_disp=0)
## get parameters back onto original scale
pp <- invlink_trans(restore(g1$par, opt_pars, fp))
print(pp)
bd <- ldmy(c("23-Mar-2020","30-Mar-2020"))
##g1R$par[4] <- -8  ## HACK/test
ed <- "1-May-2020"
r <- forecast_sim(g1$par, opt_pars,
                  ## fixed_pars = fp,
                  base_params=params,
                  start_date=min(dd$date)-15,
                  end_date=ed,
                  break_dates=bd,
                  aggregate_args = agg_list)

## aggregate_args = agg_list)
## FIXME: r can't use plot.pansim method ATM
print(ggplot(r,aes(date,value,colour=var))
      +geom_line()
      + scale_y_log10(limits=c(1,NA),oob=scales::squish)
      + geom_point(data=dplyr::mutate_at(dd,"var",trans_state_vars))
      + geom_vline(xintercept=bd,lty=2)
      )

## parameter ensemble
set.seed(101)
## HACK: NB fit is unhappy because there is severe underdispersion (because of
##  overfitting to time series); NB disp parameter is >>> 1
##  we seem to be able to get away with ignoring it completely here
##  (not needed for forecast ...)
## might help to fix it ...
e_pars <- as.data.frame(MASS::mvrnorm(200,
                                      mu=g1$par[1:4],
                                      Sigma=solve(g1$hessian[1:4,1:4])))
## tried with purrr::pmap but too much of a headache
t1 <- system.time(e_res <- plyr::alply(as.matrix(e_pars)
                                     , .margins=1
                                     , .fun=forecast_sim
                                     , aggregate_args = agg_list
                                     , base_params=params
                                     , start_date=min(dd$date)-15,
                                     , end_date=ed,
                                     , return_val="vals_only"
                                     , break_dates=bd
                                     , opt_pars = opt_pars
                                     ## , fixed_pars = fp
## breaks with fixed_pars *and* aggregate_args but OK with either?
))

## get quantiles by observation
e_res2 <- (e_res %>% bind_cols()
    %>% apply(1,quantile,c(0.1,0.5,0.9),na.rm=TRUE)
    %>% t()
    %>% as_tibble()
    %>% setNames(c("lwr","value","upr"))
)

## date/var values
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



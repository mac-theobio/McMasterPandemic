## example


library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(anytime)

first_date <- as.Date("2020-03-15")
end_date <- as.Date("2020-03-25")
dd <- data.frame(date = as.Date(first_date:end_date)
    , value = c(100,105,120,160,200,300,70,60,45,55,50)
    # , value = c(30,300,40,50,60,65,70,60,45,55,50)
    , var = "newConfirmations"
)


## 
params <- fix_pars(read_params("ICU1.csv")
                   , target=c(Gbar=6)
                   , pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
)
params[["N"]] <- 19.5e6  ## reset pop to Ontario

## breakpoints
schoolClose <- "20-Mar-2020"

bd <- anydate(c(schoolClose))
## print(bd)

opt_pars <- list(
    ## these params go to run_sim
    params=c(log_E0=4
             , log_beta0=-1
             ## fraction of mild (non-hosp) cases
             , log_mu=log(params[["mu"]])
             ## fraction of incidence reported
             ## logit_c_prop=qlogis(params[["c_prop"]]),
             ## fraction of hosp to acute (non-ICU)
             , logit_phi1=qlogis(params[["phi1"]])
             ## fraction of ICU cases dying
             ## logit_phi2=qlogis(params[["phi2"]])
    ),
    log_rel_beta0 = rep(-1, length(bd)),
    log_nb_disp=0)

t1 <- system.time(g1 <- calibrate(data=dd
                                  , base_params=params
                                  , optim_args=list(control=list(maxit=10000),hessian=TRUE)
                                  , opt_pars = opt_pars,
                                  , break_dates = bd
                                  ## , debug=TRUE
)
) ## system.time

f_args <-attr(g1,"forecast_args")
i1 <- invlink_trans(restore(g1$par,f_args$opt_pars))
params_fitted <- update(params,i1$params)
## hack to get names of params used to adjust Gbar
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

get_type <- . %>%  mutate(vtype=ifelse(var %in% c("incidence","report","newTests"),
                                       "inc","prev"))
mk_fc <- function(new_end=NULL) {
    if (!is.null(new_end)) {
        f_args$end_date <- new_end
    }
    fc <- (do.call(forecast_sim,
                   c(list(p=g1$par), f_args))
           %>% dplyr::filter(var %in% keep_vars)
           %>% get_type()
    )
    return(fc)
}
forecast1 <- mk_fc()

aa <- forecast1 %>% filter(var == "report")
(ggplot(aa,aes(x=date,y=value))
    + geom_line()
    + geom_vline(xintercept = as.Date("2020-03-20"))
    + scale_x_date(breaks="2 day")
)

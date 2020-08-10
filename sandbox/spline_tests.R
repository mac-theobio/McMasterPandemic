library(McMasterPandemic)
library(dplyr)
devtools::load_all() ## update code if necessary

## process built-in data and parameters to get something we can use
dd <- ont_all %>% trans_state_vars() %>%
       filter(var %in% c("H","report"))
params <- read_params("ICU1.csv")

## run calibrate_comb() with interesting spline settings to
##  return the time-varying beta log-linear model matrix:
X <- calibrate_comb(data=dd, params=params,
                    use_spline=TRUE,
                    spline_type="ns",
                    spline_setback=14,
                    spline_extrap="constant",
                    return="X")

unique(dd$date) # the same as rowdim of X

par(las=1)
matplot(X, type="l",lty=1, lwd=2,
        col=palette(),
        xlab="day",ylab="basis function value")


X2 <- calibrate_comb(data=dd, params=params,
                     use_spline=TRUE,
                     spline_type="ns",
                     spline_setback=14,
                     spline_extrap="constant",
                     return="args")


## same(ish), but return formula instead
ff <- calibrate_comb(data=dd, params=params,
                     use_spline=TRUE,
                     spline_type="ns",
                     return="formula")
## environment of the formula also has the model matrix in it
summary(environment(ff)$X_dat$t_vec)

ddr <- (filter(dd, var=="report")
    %>% tidyr::drop_na(value)
    %>% mutate(t_vec=as.numeric(date-min(date)))
)

## transform one-sided formula to two-sided ...
ff2 <- ff
ff2[[3]] <- ff2[[2]]
ff2[[2]] <- quote(value)
lm(ff2, data=ddr)

m <- lm(ff2, data=ddr)  ## coefficients are
pred <- X %*% coef(m)
plot(pred)

bb <- coef(m)/250  # not sure why /250?
pred2 <- X %*% bb
plot(pred2)

date <- environment(ff)$X_dat$date
stopifnot(nrow(X)==length(date))

sims <- run_sim_loglin(params=params,
               extra_pars=list(time_beta=bb),
               time_args=list(X_date=date, X=X),
               sim_args=list(start_date=date[1],end_date=tail(date,1))
               )

# run_sim_loglin(params=params[1:24],
#                extra_pars=list(time_beta=bb),
#                time_args=list(X_date=date, X=X),
#                sim_args=list(start_date=date[1],end_date=tail(date,1))
# )


head(sims)
# this gives time series of beta0
beta0_info <- attributes(sims)$params_timevar
R0 <- get_R0(params)
# relative transmission rate (which scales Rt)
beta_rel <- R0 * beta0_info$Relative_value
plot(beta_rel)

plot(sims$R)

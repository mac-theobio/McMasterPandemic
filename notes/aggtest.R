
library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(dplyr)

dd <- dplyr::filter(ont_recent,var=="newConfirmations")

dd <- (dd 
  %>% mutate(week = format(date, "%Y-%U"))
  %>% group_by(week)
  %>% summarise(date = max(date)
		, value = sum(value)
		, var = first(var)
	)
)
agg_list <- list(t_agg_start="07-Mar-2020",t_agg_period="7 days",t_agg_fun=sum)

print(ggplot(dd,aes(date,value)) + geom_point() + scale_y_log10())
## adjust parameters to sensible generation interval
params <- fix_pars(read_params("ICU1.csv"), target=c(Gbar=6),u_interval=c(-1,1),
					  pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a")))
summary(params)
g1 <- get_break_gen(data=dd, base_params=params, debug=FALSE,
						optim_args=list(control=list(maxit=10000),hessian=TRUE),
						var=if (!use_hosp) "report" else "H",
						aggregate_args = agg_list,
						debug_plot=FALSE)
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
					 end_date="1-Jun-2020",
					 break_dates=bd,
					 aggregate_args = agg_list)
## FIXME: r can't use plot.pansim method ATM
print(ggplot(r,aes(date,value,colour=var))
	 +geom_line()
	 + scale_y_log10()
	 + geom_point(data=dplyr::mutate_at(dd,"var",trans_state_vars))
	 + geom_vline(xintercept=bd,lty=2)
)

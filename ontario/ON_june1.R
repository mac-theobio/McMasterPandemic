## If you want to replicate the old pipeline, you have to go back in time
## use GH and checkout a commit on June 1st
## e58ee12f2fdae6092a660cae8201021e688a2e92




library(McMasterPandemic)
library(tidyverse)

## To make sure you are matching very thing perfectly
## load the old object and take the pieces you need to calibrate

load("inst/testdata/ONcalib_2020Jun01.rda")
ont_all_sub <- ont_cal1$mle2@data$data
params <- ont_cal1$forecast_args$base_params
# bd <- ont_cal1$forecast_args$time_args$break_dates
bd <- c("2020-03-17","2020-03-23","2020-03-28")
	
	
## This is copied over from the old script
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
					 ## fraction of ICU cases dying
					 ## logit_phi2=qlogis(params[["phi2"]])
	),
	## changes in beta at breakpoints
	logit_rel_beta0 = rep(-1, length(bd)),
	## NB dispersion
	log_nb_disp=NULL)

## do the calibration
t_ont_cal1 <- system.time(new_fit <- calibrate(data=ont_all_sub
																								, base_params=params
																								, opt_pars = opt_pars
																								, time_args=list(break_dates = bd)
)
) 
##### matt's code

## Use my new calibrated object instead
# load("refactor/ontario-rda/ONcalib_2021May17.rda")
end_date <- as.Date(new_fit$mle2@data$end_date) + 30
prediction_new <- predict(new_fit,end_date=end_date) %>% filter(! (var %in% c('cumRep')))

load("inst/testdata/ONcalib_2020Jun01.rda")
end_date <- as.Date(ont_cal1$mle2@data$end_date) + 30
prediction_old <- predict(ont_cal1,end_date=end_date) %>% filter(! (var %in% c('cumRep')))

replace_na <- function(x) {x[is.na(x)]<- 0; x }
stopifnot(all(prediction_new$date == prediction_old$date))
stopifnot(all(prediction_new$var == prediction_old$var))
stopifnot(all(prediction_new$vtype == prediction_old$vtype))
merged <- dplyr::inner_join(prediction_new, prediction_old, by=c("date", "var", "vtype"), suffix=c("_new", "_old")) %>%
	mutate(value_new = replace_na(value_new), value_old = replace_na(value_old)) %>%
	mutate(abs_diff = abs(value_new-value_old))


ggplot_df <- merged %>% filter(var=='incidence')
ggplot(ggplot_df) +
	theme(text=element_text(size=20))+
	geom_line(aes(x=date, y=value_new, color='new'), size=1) +
	geom_line(aes(x=date, y=value_old, color='old'), size=1) +
	labs(x='date', y='incidence', title='New vs old prediction discrepancies', col='')


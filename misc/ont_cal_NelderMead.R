
library(McMasterPandemic)
library(tidyverse)

dir.create(file.path('weekly_calibrations'), showWarnings = FALSE)

keep_vars <- c("H", "ICU", "death", "report")
## data since 15 March

ont_all_sub <- (ont_all
                %>% mutate_at("var", trans_state_vars)
                %>% filter(var %in% keep_vars)
)

## adjust mean GI
params <- fix_pars(read_params("ICU1.csv"),
                   target = c(Gbar = 6),
                   pars_adj = list(c("sigma", "gamma_s", "gamma_m", "gamma_a"))
)
params[["N"]] <- 14.57e6 ## reset pop to Ontario



bd <- c("2020-03-17", "2020-03-23", "2020-03-28")


## This is copied over from the old script
opt_pars <- list(
  ## these params are part of the main parameter vector: go to run_sim()
  params = c(
    log_E0 = 4 ## initial exposed
    , log_beta0 = -1 ## initial baseline transmission
    ## fraction of mild (non-hosp) cases
    , log_mu = log(params[["mu"]])
    ## fraction of incidence reported
    ## logit_c_prop=qlogis(params[["c_prop"]]),
    ## fraction of hosp to acute (non-ICU)
    , logit_phi1 = qlogis(params[["phi1"]])
    ## fraction of ICU cases dying
    ## logit_phi2=qlogis(params[["phi2"]])
  ),
  ## changes in beta at breakpoints
  logit_rel_beta0 = rep(-1, length(bd)),
  ## NB dispersion
  log_nb_disp = NULL,
  time_params = c(1, 1, 1)
)


timevar_df <- data.frame(Date = bd,
                         Symbol = c("beta0", "beta0", "beta0"),
                         Relative_value = c(NA, NA, NA))

start_time <- Sys.time()
ont_cal1 <- calibrate(
  data = ont_all_sub,
  base_params = params,
  opt_pars = opt_pars,
  time_args = list(params_timevar = timevar_df),
  sim_args = list(step_args = list())
)

end_time <- Sys.time()

time_diff <- end_time - start_time
print(time_diff)

file = file.path("../weekly_calibrations", paste("ONcalib_", format(Sys.time(), "%Y%b%d"), ".rda", sep=""))
save("ont_cal1", "bd", "ont_all_sub", "time_diff", file=file)

end_date = as.Date(ont_cal1$mle2@data$end_date) + 30
prediction <- predict(ont_cal1, end_date = end_date) %>%
  filter(!(var %in% c("cumRep")))


ggplot_df <- prediction %>% filter(var == "incidence")

file = file.path("../weekly_calibrations", paste("prediction_", format(Sys.time(), "%Y%b%d"), ".png", sep=""))
ggplot(ggplot_df) +
	theme(text = element_text(size = 20)) +
	geom_line(aes(x = date, y = value), size = 1) +
	labs(x = "date", y = "incidence", title = "", col = "")
	ggsave(file=file, width = 16, height = 10)

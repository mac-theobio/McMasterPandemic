
library(McMasterPandemic)
library(tidyverse)


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
  log_nb_disp = NULL
)


timevar_df <- data.frame(Date = bd,
                         Symbol = c("beta0", "beta0", "beta0"),
                         Relative_value = c(NA, NA, NA))

no_hazard <- calibrate(
  data = ont_all_sub,
  base_params = params,
  opt_pars = opt_pars,
  time_args = list(params_timevar = timevar_df),
  sim_args = list(step_args = list(do_hazard = FALSE))
)

do_hazard <- calibrate(
  data = ont_all_sub,
  base_params = params,
  opt_pars = opt_pars,
  time_args = list(params_timevar = timevar_df),
  sim_args = list(step_args = list(do_hazard = TRUE))
)

end_date <- as.Date(no_hazard$mle2@data$end_date) + 30
prediction_no_hazard <- predict(no_hazard, end_date = end_date) %>%
  filter(!(var %in% c("cumRep")))

end_date <- as.Date(do_hazard$mle2@data$end_date) + 30
prediction_do_hazard <- predict(do_hazard, end_date = end_date) %>%
  filter(!(var %in% c("cumRep")))

replace_na <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}
stopifnot(all(prediction_do_hazard$date == prediction_no_hazard$date))
stopifnot(all(prediction_do_hazard$var == prediction_no_hazard$var))
stopifnot(all(prediction_do_hazard$vtype == prediction_no_hazard$vtype))

merged <- dplyr::inner_join(prediction_do_hazard,
                            prediction_no_hazard,
                            by = c("date", "var", "vtype"),
                            suffix = c("_do_hazard", "_no_hazard")) %>%
  mutate(value_do_hazard = replace_na(value_do_hazard),
         value_no_hazard = replace_na(value_no_hazard)) %>%
  mutate(abs_diff = abs(value_do_hazard - value_no_hazard))


ggplot_df <- merged %>% filter(var == "incidence")
ggplot(ggplot_df) +
  theme(text = element_text(size = 20)) +
  geom_line(aes(x = date, y = value_do_hazard, color = "do_hazard=TRUE"), size = 1) +
  geom_line(aes(x = date, y = value_no_hazard, color = "do_hazard=FALSE"), size = 1) +
  labs(x = "date", y = "incidence", title = "do_hazard prediction discrepancies", col = "")
ggsave("discrepancies.png", width = 16, height = 10)


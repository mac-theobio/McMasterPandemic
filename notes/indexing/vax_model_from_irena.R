library(McMasterPandemic)

start_date <- "2020-02-01"
end_date <- "2020-09-01"

## initialize params
base_params <- read_params("PHAC.csv")

## what is the structure of the base params?
str(base_params)

base_sim <- run_sim(
  params = base_params,
  start_date = start_date,
  end_date = end_date
)

plot(base_sim)

vax_params <- expand_params_vax(
  params = base_params,
  model_type = "twodose"
)

setdiff(names(vax_params), names(base_params))

base_state <- make_state(params = base_params)

attr(base_state, "epi_cat")

vax_state <- expand_state_vax(
  x = base_state,
  model_type = "twodose",
  unif = FALSE
)

vax_state

options(macpan_pfun_method = "grep")
vax_sim <- run_sim(
  params = vax_params,
  start_date = start_date,
  end_date = end_date
)

plot(vax_sim)

vax_sim_full <- run_sim(
  params = vax_params,
  start_date = start_date,
  end_date = end_date,
  condense_args = list(keep_all = TRUE)
)

names(vax_sim_full)

##' @param x a \code{pansim} object generated using the model with vaccination
plot_vaxified_sim <- function(x){
  p <- (x
        %>% select(-!(contains("_")|date))
        %>% pivot_longer(cols = -date)
        %>% tidyr::separate(
          col = name,
          into = c("variable", "vax_cat")
        )
        %>% mutate(variable = forcats::as_factor(variable),
                   vax_cat = forcats::as_factor(vax_cat)) ## to keep states and vax_cat ordered by disease history instead of alphabetically
        %>% ggplot(aes(x = date, y = value,
                       colour = variable))
        + geom_line()
        + facet_grid(rows = vars(variable),
                     cols = vars(vax_cat),
                     scales = "free_y")
        + theme(axis.title = element_blank())
  )
  return(p)
}

plot_vaxified_sim(vax_sim_full)

vax_sim_full2 <- run_sim(
  params = vax_params,
  start_date = start_date,
  end_date = end_date,
  condense_args = list(keep_all = TRUE),
  ## allocate 50% of the doses administered per day to second doses 30 days after the simulation start date
  params_timevar = data.frame(
    Date = as.Date(start_date) + 30,
    Symbol = "vax_prop_first_dose",
    Value = 0.5,
    Type = "rel_orig" # take 0.5 of the parameter in the *original* params list (which happens to be 1 for this parameter) as opposed to taking it relative to its previous value ("rel_prev") in params_timevar (there isn't one)
  )
)

params_timevar <- data.frame(
  Date = c(as.Date(start_date) + 30,
           as.Date(start_date) + 60),
  Symbol = c("vax_prop_first_dose", "beta0"),
  Value = rep(0.5, 2),
  Type = rep("rel_orig", 2)
)

## generate reports from sim
synth_reports <- (run_sim(
  params = vax_params,
  start_date = start_date,
  end_date = end_date,
  # do the same thing with the switch to second doses as above, but now also cut the original transmission rate to 50% of its value 60 days after the simulation start date
  params_timevar = params_timevar
)
## reshape into the correct format for input data passed to calibrate()
%>% mutate(value=round(report), var="report")
%>% select(date, value, var)
%>% na.omit()
)

## set up optimization parameters
## (base parameter values)
opt_pars <- list(
  params = c(beta0 = 0.6), ## set initial guess for beta0
  time_params = c(0.8) ## initial guess for change in beta0 on the one and only break date (guess 80% of the base value)
)
## (time-varying relative values: insert an NA wherever you want a parameter to be calibrated)
params_timevar_calib <- (params_timevar
                         %>% mutate(Value = ifelse(Symbol == "beta0",
                                                   NA,
                                                   Value))
)


fitted_mod <- calibrate(
  base_params = vax_params,
  data = synth_reports,
  debug = TRUE,
  opt_pars = opt_pars,
  time_args = list(
    params_timevar = params_timevar_calib
  ),
  sim_args = list(
    ndt = 1,
    step_args = list(do_hazard = TRUE)
  ) ## there are both the defaults currently but i'm putting them here to emphasize that for the vaxified model's current implementation,
  ## we need ndt = 1 because the vax rate is specified as doses *per day*, implicitly assuming that the simulation algorithm takes daily steps
  ## and we need do_hazard = TRUE to avoid accidentally stepping into negative state space when vaccination occurs relatively quickly
)





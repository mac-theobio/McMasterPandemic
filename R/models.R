##' Represent a Standard Model as a flexmodel
##'
##' @export
make_unflexmodel <- function(params,
                             state = NULL,
                             start_date = "2020-03-20",
                             end_date = "2020-05-1",
                             params_timevar = NULL,
                             step_args = list()) {

  ## check_start = Sys.time()
  spec_check("0.0.5", "run_sim with TMB")

  ## may need to modify this when we start updating time-varying parameters
  ## on the c++ side (currently params0 should always equal params)
  params0 <- params
  state0 <- state

  step_args_to_flex <- c("do_hazard")
  flex_args <- c(
    list(
      params = params,
      state = state,
      start_date = start_date,
      end_date = end_date,
      params_timevar = params_timevar
    ),
    step_args[names(step_args) %in% step_args_to_flex]
  )
  init_end <- Sys.time()

  ## make_model_start = Sys.time()
  model <- (init_model
            %>% do.call(flex_args)
            %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
            %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
            %>% add_rate("Ia", "R", ~ (gamma_a))
            %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
            %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
            %>% add_rate("Im", "R", ~ (gamma_m))
            %>% add_rate("Is", "H", ~
                           (1 - nonhosp_mort) * (phi1) * (gamma_s))
            %>% add_rate("Is", "ICUs", ~
                           (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
            %>% add_rate("Is", "ICUd", ~
                           (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
            %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
            %>% add_rate("ICUs", "H2", ~ (psi1))
            %>% add_rate("ICUd", "D", ~ (psi2))
            %>% add_rate("H2", "R", ~ (psi3))
            %>% add_rate("H", "R", ~ (rho))
            %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
            %>% add_rate("S", "E", ~
                           (Ia) * (beta0) * (1 / N) * (Ca) +
                           (Ip) * (beta0) * (1 / N) * (Cp) +
                           (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
                           (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
            %>% add_parallel_accumulators(c("X", "N", "P", "V"))
            %>% add_tmb_indices()
  )
  return(model)
}

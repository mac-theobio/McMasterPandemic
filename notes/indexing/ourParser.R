# Use our own "parser" proposed in tmb_data_struct.md

library(TMB)
compile("ourParser.cpp")
dyn.load(dynlib("ourParser"))
#library("McMasterPandemic")

source("utilities.R")

rs = mk_ratemat_struct(
  rate("E", "Ia", ~ (alpha) * (sigma)),
  rate("E", "Ip", ~ (1 - alpha) * (sigma)),
  rate("Ia", "R", ~ (gamma_a)),
  rate("Ip", "Im", ~ (mu) * (gamma_p)),
  rate("Ip", "Is", ~ (1 - mu) * (gamma_p)),
  rate("Im", "R", ~ (gamma_m)),
  rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s)),
  rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s)),
  rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s)),
  rate("Is", "D", ~ (nonhosp_mort) * (gamma_s)),
  rate("ICUs", "H2", ~ (psi1)), ## ICU to post-ICU acute care
  rate("ICUd", "D", ~ (psi2)), ## ICU to death
  rate("H2", "R", ~ (psi3)), ## post-ICU to discharge
  ## H now means 'acute care' only; all H survive & are discharged
  #list(from = "H", to = "D", formula = ~ 0),
  rate("H", "R", ~ (rho)), ## all acute-care survive
  rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s)), ## assuming that hosp admissions mean *all* (acute-care + ICU)
  # force of infection
  rate("S",  "E", ~
         (Ia) * (beta0) * (1/N) * (Ca) +
         (Ip) * (beta0) * (1/N) * (Cp) +
         (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
         (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
)

ds = to_tmb(rs)

print(ds)

# Ready to go
dd <- MakeADFun(data = list(state = c(state),
                            ratemat = as(M, "dgTMatrix"),
                            from = ds$from,
                            to = ds$to,
                            count = ds$count,
                            spi = ds$spi,
                            modifier = ds$modifier,
                            dt = 1,
                            do_hazard = TRUE,
                            stoch_proc = FALSE,
                            do_exponential = FALSE,
                            testwt_scale = "N"
                            ),
                parameters = list(params=c(params)),
                DLL='ourParser'
               )

new_state = dd$report()$state

new_state



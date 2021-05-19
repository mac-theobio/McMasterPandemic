devtools::load_all()
library(dplyr)
library(ggplot2)
library(patchwork)

## sim settings
start_date <- "2020-02-01"
end_date <- "2020-06-01"

time_pars <- data.frame(Date=c(start_date, "2020-03-01"),
                        Symbol="vax_doses_per_day",
                        Relative_value = c(0, 1e5))

## params
## set up vax params
vax_doses_per_day <- 1

## initialize params
base_params <- update(read_params("PHAC.csv")
                      , N = 1.4e7
)

vax_params <- update(expand_params_vax(base_params,
                                       vax_doses_per_day = vax_doses_per_day))

## initialize states
vax_state <- make_state(params = vax_params)

res_vax <- run_sim(vax_params, vax_state,
                   start_date = start_date, end_date = end_date,
                   params_timevar = time_pars,
                   step_args = list(do_hazard = TRUE),
                   condense_args = list(keep_all = TRUE))

p1 <- ggplot(res_vax, aes(x = date, y = report)) + geom_point()

p2 <- (ggplot(get_doses_per_day(res_vax),
              aes(x = date, y = total_doses_per_day))
       + geom_point()
       + geom_vline(aes(xintercept = as.Date(time_pars$Date[2])),
                    linetype = "dashed",
                    colour = "grey60")
       + annotate("text",
                  x = as.Date(time_pars$Date[2])-8,
                  y = 0.5*time_pars$Relative_value[2],
                  angle = 90,
                  colour = "grey60",
                  label = paste0("vaccination starts\nat rate of\n",
                                 format(time_pars$Relative_value[2], scientific = FALSE, big.mark = ","), " doses/day"))
)

ggsave(filename = "git_push/demo_vaccination_switch.png",
       plot = (p1 / p2))

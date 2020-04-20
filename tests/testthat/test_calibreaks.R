library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)

params <- fix_pars(read_params("ICU1.csv"))


start_date <- anydate("2020-01-01")
end_date <- anydate("2020-04-01")
break1 <- "2020-02-22"
bd <- anydate(c(break1))
rel_break1 <- 0.3

params[["N"]] <- 1e7

sim1S <- run_sim_break(params
   , start_date=start_date
   , end_date=end_date
   , break_dates = bd
   , rel_beta0= rel_break1
)

simI <- (sim1S
   %>% transmute(date
      , I = S*foi 
      )         
)

## Is it suppose to look like this?

print(ggplot(simI,aes(x=date,y=I))
   + geom_line()
   + geom_vline(xintercept = bd)
)

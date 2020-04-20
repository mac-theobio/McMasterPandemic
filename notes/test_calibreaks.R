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
plot(sim1S)

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

pp <- pivot(sim1S) %>% filter(! var %in% c("S","t"))
print(gg_all <- ggplot(pp,aes(x=date,y=value,colour=var)) 
   + geom_line()
   + geom_vline(xintercept = bd)
   + scale_y_log10()
)


sim1S2 <- run_sim_break(params
   , start_date=start_date
   , end_date=end_date
   , break_dates = bd
   , rel_beta0= rel_break1
    , ndt=20
     )

pp2 <- (pivot(sim1S2)
    %>% filter(! var %in% c("S","t"))
    %>% mutate(foi=(var=="foi"))
    %>% filter(date >= as.Date("2020-02-20"), date <= as.Date("2020-03-01"))
)
print(gg_all  %+% pp2
      + geom_point()
      + facet_wrap(~foi,scale="free_y")
)

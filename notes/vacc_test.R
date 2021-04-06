library(McMasterPandemic)
library(cowplot)
library(shellpipes)
library(ggplot2)
theme_set(theme_bw())

commandEnvironments()
makeGraphics()


devtools::load_all("..")
p0 <- fix_pars(read_params("ICU1.csv"), target = c(R0=1.5,Gbar=6))
p1 <- update(p0, vacc=0.001)
p2 <- update(p0, vacc=1e-7)

r1 <- run_sim(params=p0,
              start_date="01-Mar-2020",
              end_date="01-Jul-2020")

r2 <- run_sim(params=p1,
              start_date="01-Mar-2020",
              end_date="01-Jul-2020")

tv <- data.frame(Date=c("15-Mar-2020"
    , "30-Mar-2020"
    )
    , Symbol=rep("vacc",1)
    , Relative_value=c(1e5,2e5)
)
r3 <- run_sim(params=p2,
              start_date="01-Mar-2020",
              end_date="01-Jul-2020",
              params_timevar=tv)
r4 <- run_sim(params=update(p2,vacc=1e-2),
              start_date="01-Mar-2020",
              end_date="01-Jul-2020")


print(plot_grid(plot(r1,keep_states=c("S","R"))
	, plot(r2,keep_states=c("S","R")) 
	, plot(r3,keep_states=c("S","R")) 
	, plot(r4,keep_states=c("S","R")) 
	, ncol=2)
)

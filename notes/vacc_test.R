library(McMasterPandemic)
library(cowplot)
devtools::load_all("..")
p0 <- read_params("ICU1.csv")
p1 <- update(p0, vacc=0.001)
p2 <- update(p0, vacc=1e-7)

r1 <- run_sim(params=p0,
              start_date="01-Mar-2020",
              end_date="01-Jul-2020")

r2 <- run_sim(params=p1,
              start_date="01-Mar-2020",
              end_date="01-Jul-2020")

tv <- data.frame(Date=c("01-Apr-2020",
                        "30-Apr-2020"),
                 Symbol=rep("vacc",2),
                 Relative_value=c(1e5,1e6))
r3 <- run_sim(params=p2,
              start_date="01-Mar-2020",
              end_date="01-Jul-2020",
              params_timevar=tv)
r4 <- run_sim(params=update(p2,vacc=1e-7),
              start_date="01-Mar-2020",
              end_date="01-Jul-2020")


plot_grid(
    plot(r1,keep_states=c("S","R")),
    plot(r2,keep_states=c("S","R")),
    plot(r3,keep_states=c("S","R")),
    plot(r4,keep_states=c("S","R")),
    ncol=2)

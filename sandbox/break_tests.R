## spikes

params <- read_params("ICU1.csv")
r1 <- run_sim_break(params, time_args=list(break_dates="2020-03-01"),
                    sim_args=list(start_date="2020-02-01", end_date="2020-04-01"),
                    extra_pars=list(rel_beta0 = 0.5))

r2 <- run_sim_break2(params, time_args=list(break_dates="2020-03-01"),
                    sim_args=list(start_date="2020-02-01", end_date="2020-04-01"),
                    extra_pars=list(rel_beta0 = 0.5))



## changing other parameters

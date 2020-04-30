library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(anytime)

print(coef(g1))
plot(run_sim(
	params=update(coef(g1), c(obs_disp=20, proc_disp=1))
	, stoch=c(proc=TRUE,obs=FALSE)
	, stoch_start = c(proc=anydate("2020-04-12"), obs=anydate("2020-04-02"))
), log=TRUE)

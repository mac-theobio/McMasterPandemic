library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(anytime)

options(stringsAsFactors=FALSE)
res <- readRDS("ont_fac_1011_fit.rds")
dat <- res$data

set.seed(2121)
end <- "2020-06-20"
nsim <- 50
epistart <- min(dat$date) - 15
data_end <- max(dat$date)
print(data_end) 

## Plot the fit and look at error assumptions


somevars <- predict(object=res$fit, ensemble=TRUE
   , end_date=end
   , nsim=nsim
#  , stoch=c(proc=TRUE,obs=TRUE)
#  , new_params=c(obs_disp=20, proc_disp=1.0)
#  , stoch_start = c(proc=data_end, obs=epistart)
)
plot(res$fit, data=dat)             
plot(somevars, data=dat)

library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(anytime)

options(stringsAsFactors=FALSE)

res_list <- readRDS("fit_cali.rds")

cali_mod <- res_list[[1]]

dat <- res_list[[2]]

set.seed(2121)
end <- "2020-06-20"
nsim <- 50
epistart <- min(dat$date) - 15
data_end <- max(dat$date)
print(data_end) 

## Plot the fit and look at error assumptions

## devtools::load_all("..")
somevars <- predict(object=cali_mod, ensemble=TRUE
        , stoch=c(proc=TRUE,obs=TRUE)
	, new_params=c(obs_disp=20 , proc_disp=1.0)
	, end_date=end
	, stoch_start = c(proc=data_end+1, obs=epistart)
	, nsim=nsim
          )
##        , keep_vars="Rt"

plot(cali_mod, data=dat)					
plot(somevars, data=dat)

## Rt calculation is enabled by the presence of "Rt" in keep_vars
Rtenv <- predict(object=cali_mod, ensemble=TRUE
               , keep_vars="Rt"                 
               , nsim=nsim
                 )
## compare against free-standing version
plot(Rtenv) + geom_line(data=get_Rt(cali_mod),aes(y=R0))


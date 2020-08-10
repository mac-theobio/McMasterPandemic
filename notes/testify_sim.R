library(tidyverse)
source("makestuff/makeRfuns.R")
## commandEnvironments()

## WHICH of these do you want today?
## library(devtools); load_all("../")
library("McMasterPandemic")

## if (interactive()) {
use_ode <- TRUE
testwt_scale <- "none" ## or "N" or "sum_u"
## }
## Magic at the beginning
set.seed(0807)
start <- as.Date("2020-01-01")
end <- as.Date("2020-06-01")

fn <- if (interactive()) "PHAC_testify.csv" else matchFile(".csv$")
params <- (read_params(fn)
    %>% fix_pars(target=c(R0=2.5, Gbar=6))
    %>% update(
            N=1.5e7       ## population of Ontario
        )
)

paramsw0 <- params[!grepl("^W",names(params))] ## removing all of the regular W-parameters
class(paramsw0) <- "params_pansim"

print(paramsw0)
summary(paramsw0)

testing_intensity <- c(0.002, 0.02, 0.2)
W_asymp <- c(0.01, 0.1,1)
iso_t <- c(0,0.5,0.9,1)

simlist <- list()
for(i in W_asymp) {
    for (j in iso_t) {
        for (k in testing_intensity) {
            cat(i,j,k,"\n")
            paramsw0 <- update(paramsw0, W_asymp=i, iso_t = j, testing_intensity=k)
            sims <- (run_sim(params = paramsw0, ratemat_args = list(testify=TRUE)
                           , start_date = start
                           , end_date = end
                           , use_ode = use_ode
                           , step_args = list(testwt_scale=testwt_scale)
                             ##			, condense_args=list(keep_all=TRUE) checkout the expanded version
                             )
                %>% mutate(W_asymp = i
                         , iso_t = j
                         , testing_intensity=k
                           )
            )
            simlist <- c(simlist,list(sims))
        } ## loop over testing intensity
    } ## loop over iso_t
} ## loop over W_asymp
simframe <- bind_rows(simlist)

print(simframe)

simdat <- (simframe
    %>% transmute(date
                , incidence
                , postest
                , total_test = postest + negtest
                , positivity = postest/total_test
                , report
                , W_asymp
                , iso_t
                , testing_intensity
                  )
    %>% gather(key="var",value="value",-c(date, W_asymp, iso_t, testing_intensity))
)

saveVars(simdat, params)

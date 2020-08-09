library(tidyverse)
source("makestuff/makeRfuns.R")
## commandEnvironments()

## WHICH of these do you want today?
## library(devtools); load_all("../")
library("McMasterPandemic")

## Magic at the beginning
set.seed(0807)
start <- as.Date("2020-01-01")
end <- as.Date("2020-06-01")

fn <- if (interactive()) "PHAC_testify.csv" else matchFile(".csv$")
pp <- (read_params(fn)
    %>% update(
            iso_t=0         ## no isolation of tested
            , N=1.5e7       ## population of Ontario
            , beta0=0.4     ## Make a semi-realistic R0
            , testing_intensity=0.2  ## high testing intensity so we can see what's going on
        )
)

ppw0 <- pp[!grepl("^W",names(pp))] ## Copying BMB, removing all of the regular W-parameters
class(ppw0) <- "params_pansim"

print(ppw0)
summary(ppw0)

W_asymp <- c(0.01, 0.1,1)
iso_t <- c(0,0.5,0.9,1)

simlist <- list()
for(i in W_asymp) {
    for(j in iso_t) {
        ppw0 <- update(ppw0, W_asymp=i, iso_t = j)
        sims <- (run_sim(params = ppw0, ratemat_args = list(testify=TRUE)
                       , start_date = start
                       , end_date = end
                         ##			, condense_args=list(keep_all=TRUE) checkout the expanded version
                         )
            %>% mutate(W_asymp = i,
                       iso_t = j)
        )
        simlist <- c(simlist,list(sims))
    }
}
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
                  )
    %>% gather(key="var",value="value",-c(date, W_asymp, iso_t))
)

saveVars(simdat)

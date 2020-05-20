##
## possible params to set/command-line arguments:
##  Gbar = 6
##  which vars to fit?
library(McMasterPandemic)
library(tidyverse)

## clean Ontario data are automatically loaded
## load("ontario_clean.RData")
## translate variable names to internally used values
## drop unused variables
keep_vars <- c("H","ICU","death","report")
## data since 15 March
ont_recent_sub <- (ont_recent
    %>% mutate_at("var",trans_state_vars)
    %>% filter(var %in% keep_vars)
)


ont_all_sub <- (ont_all
    %>% mutate_at("var",trans_state_vars)
    %>% filter(var %in% keep_vars)
)

## unique(ont_recent_sub$var)

## adjust mean GI
params <- fix_pars(read_params("ICU1.csv")
    , target=c(Gbar=6)
    , pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
)
params[["N"]] <- 14.57e6  ## reset pop to Ontario


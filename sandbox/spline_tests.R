library(McMasterPandemic)
library(dplyr)
devtools::load_all() ## update code if necessary
dd <- ont_all %>% trans_state_vars() %>%
       filter(var %in% c("H","report"))
params <- read_params("ICU1.csv")

## 
X <- calibrate_comb(data=dd, params=params,
                    use_spline=TRUE,
                    spline_type="ns",
                    spline_setback=14,
                    spline_extrap="constant",
                    return="X")
matplot(X, ylab="")

ff <- calibrate_comb(data=dd, params=params,
                     use_spline=TRUE,
                     spline_type="ns",
                     return="formula")
## environment is already loaded
summary(environment(ff)$X_dat$t_vec)



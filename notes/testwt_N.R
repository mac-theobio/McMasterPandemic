library(shellpipes)

commandEnvironments()

use_ode <- FALSE
testwt_scale <- "N"  ## experiment with rescaling testing weights
testing_type <- c("logistic") ## we know constant is good, linear is kind of like logistic
testing_type <- c("constant")
testing_time <- "sample"
iso_t <- c(0,1)
omega <- c(1)
start <- as.Date("2020-01-01")
end <- as.Date("2021-01-01")
pop <- 1.5e7
R0 <- c(1.5,2)
W_asymp <- 0.05
Gbar <- c(6)
set.seed(1104)
keep_vars <- c("death/H/postest")
keep_all <- FALSE
stoch_obs <- FALSE
min_testing <- 0.002
max_testing <- 0.01
opt_testify <- FALSE


testing_intensity <- min_testing
testing_intensity <- c(0.002,0.01)

saveEnvironment()


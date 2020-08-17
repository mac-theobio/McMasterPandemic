use_ode <- FALSE
testwt_scale <- "N"  ## experiment with rescaling testing weights
# testing_intensity <- c(0.002, 0.005, 0.01)
W_asymp <- c(0.01, 0.1,1)
iso_t <- c(0,0.5,0.9,1)
start <- as.Date("2020-01-01")
end <- as.Date("2020-10-01")
pop <- 1.5e7
R0 <- 2.5
Gbar <- 6
set.seed(0807)
keep_vars <- c("hosp/death/report")
constant_testing <- c(TRUE, FALSE)
keep_all <- FALSE

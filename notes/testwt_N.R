use_ode <- FALSE
testwt_scale <- "N"  ## experiment with rescaling testing weights
testing_type <- c("constant","linear","logistic") ## we know constant is good, linear is kind of like logistic
testing_type <- c("linear")
testing_intensity <- c(0.002)
iso_t <- c(0,0.5,0.9,1)
omega <- c(0.2,1)
start <- as.Date("2020-01-01")
end <- as.Date("2020-10-01")
pop <- 1.5e7
R0 <- 2.5
Gbar <- 6
set.seed(0807)
keep_vars <- c("death/H/postest")
keep_all <- FALSE

min_testing <- 0.001 ##3e-5
max_testing <- 0.01 ##2e-3


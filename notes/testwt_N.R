use_ode <- FALSE
testwt_scale <- "N"  ## experiment with rescaling testing weights
testing_intensity_type <- c("constant","linear","logistic")
testing_intensity <- c(0.002)
iso_t <- c(0,0.5,0.9,1)
omega <- c(0.2,1)
start <- as.Date("2020-03-01")
end <- as.Date("2020-08-01")
pop <- 1.5e7
R0 <- 2.5
Gbar <- 6
set.seed(0807)
keep_vars <- c("hosp/death/report")
constant_testing <- c(TRUE, FALSE)
keep_all <- FALSE

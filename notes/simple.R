use_ode <- FALSE
testwt_scale <- "N" ## or "N" or "sum_u"
testing_intensity <- c(0.01)
testing_time="sample"
W_asymp <- c(0.1, 1)
iso_t <- c(0, 1)
start <- as.Date("2020-01-01")
end <- as.Date("2020-12-01")
pop <- 1.5e7
R0 <- 2
Gbar <- 12
set.seed(0807)
keep_all <- FALSE


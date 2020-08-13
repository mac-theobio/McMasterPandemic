use_ode <- FALSE
testwt_scale <- "none" ## or "N" or "sum_u"
testing_intensity <- c(0.002, 0.02, 0.2)
W_asymp <- c(0.01, 0.1,1)
iso_t <- c(0,0.5,0.9,1)
start <- as.Date("2020-01-01")
end <- as.Date("2020-10-01")
pop <- 1.5e7
R0 <- 2.5
Gbar <- 6
set.seed(0807)
keep_all <- FALSE

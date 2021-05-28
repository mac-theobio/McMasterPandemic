library(McMasterPandemic)
library(dplyr)
testLevel <- if (nzchar(s <- Sys.getenv("MACPAN_TEST_LEVEL"))) as.numeric(s) else 1

## sim example
params <- read_params("ICU1.csv")
paramsS <- update(params, c(proc_disp = 0.1, obs_disp = 100))
state <- make_state(params = params)
sdate <- "2020-02-10" ## arbitrary!
set.seed(101)
res1 <- run_sim(params, state, start_date = sdate, end_date = "2020-06-01")
res1_S <- update(res1, params = paramsS, stoch = c(obs = TRUE, proc = TRUE))

cdat <- (res1_S
    %>% pivot()
    %>% filter(
        var == "report",
        date > as.Date("2020-03-01"),
        date < as.Date("2020-04-15")
    )
)

pt <- data.frame(
    Date = c("2020-03-23", "2020-03-30", "2020-04-01"),
    Symbol = rep("beta0", 3), Relative_value = c(0, NA, 0)
)

vague_priors <- list(~ dlnorm(time_params[1], meanlog = -1, sd = 5))
strong_priors <- list(~ dlnorm(time_params[1], meanlog = -1, sd = 0.01))

curve(dlnorm(x, meanlog = -1, sd = 0.5), from = 0, to = 1.5)
abline(v = plogis(1))

cfun <- function(priors = NULL, debug_plot = FALSE, debug = FALSE) {
    calibrate(
        data = cdat, base_params = params,
        time_args = list(params_timevar = pt),
        opt_pars = list(
            params = c(log_E0 = 4, log_beta0 = -1),
            log_time_params = -1,
            log_nb_disp = NULL
        ),
        priors = priors,
        debug_plot = debug_plot,
        debug = debug
    )
}

if (testLevel > 1) {
    c0 <- cfun()
    c1 <- cfun(priors = vague_priors)
    c2 <- cfun(priors = strong_priors)

    cList <- list(no_prior = c0, vague_priors = c1, strong_priors = c2)
    param_tab <- t(sapply(cList, function(x) unlist(coef(x, "fitted"))))

    print(param_tab)
}

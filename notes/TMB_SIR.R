library(TMB)
library(bbmle)
compile("SIR.cpp")
dyn.load(dynlib("SIR"))
library(McMasterPandemic)

very_fake_data = ceiling(dnorm(seq(from = -3, to = 3, length.out = 100)) * 100)

sim_obj = MakeADFun(
  data = list(
    obs_incidence = very_fake_data,  # only required to create simulation object
    N = 1e7,  # total number of people
    beta = as(Matrix::sparseMatrix(i = 1, j = 1, x = 0), "dgTMatrix") # trying to break AD
  ),
  parameters = list(
    log_beta = log(0.05),
    log_gamma = log(0.01),
    #log_I0 = log(0.06),
    log_nbdisp = log(20)
  ),
  DLL = "SIR"
)

r_nan = function(x, v) {x[is.nan(x)] = v; return(x)}
d_nan = function(x, v) x[!is.nan(x)]

set.seed(1)
sims = sim_obj$simulate()
simulated_incidence = d_nan(sims$obs_incidence, 0)
underlying_incidence = d_nan(exp(sims$log_incidence), 0)
plot(underlying_incidence,
     type = 'l', las = 1,
     ylim = range(c(underlying_incidence, simulated_incidence)))
points(simulated_incidence)

set.seed(1)
obj = MakeADFun(
  data = list(
    obs_incidence = simulated_incidence,
    N = 1e7,
    beta = as(Matrix::sparseMatrix(i = 1, j = 1, x = 0), "dgTMatrix")
  ),
  parameters = list(
    log_beta = log(0.05) + rnorm(1, sd = 10),
    log_gamma = log(0.0007) + rnorm(1, sd = 10),
    #log_I0 = log(0.06) + rnorm(1, sd = 0.1),
    log_nbdisp = log(20) + rnorm(1, sd = 10)
  ),
  DLL = "SIR"
)

obj$fn()
obj$gr()
obj$he()

obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt$convergence
opt$par
sim_obj$par

obj$fn(opt$par)
obj$gr(opt$par)
obj$he(opt$par)

cov2cor(solve(opt$hessian))

fn = function(log_beta, log_gamma, log_I0, log_nbdisp) {
  obj$fn(log_beta, log_gamma, log_I0, log_nbdisp)
}
mm = mle2(fn, obj$par)
mm
sim_obj$par

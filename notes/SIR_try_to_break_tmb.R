library(TMB)
compile("SIR_try_to_break_tmb.cpp")
dyn.load(dynlib("SIR_try_to_break_tmb"))

# ------------------------------------------
# Trying to break TMB Automatic
# Differentiation (AD)
# ------------------------------------------
# We introduce two ways to represent
# transmission:
#  1. y = log(beta) = log(x - 1)
#  2. x = beta + 1 = exp(y) + 1
#
# Store y in a TMB parameter
# Store x in a 1x1 sparse data matrix
#
# On the C++ side, at every function
# evaluation start by storing exp(y) + 1
# in x. Then use x in expressions
# required to compute the likelihood.
#
# This approach should fail to give
# reasonable/correct gradients if
# the AD machinery is not implemented
# for this sparse matrix class.
#
# Failure would clearly indicate
# that we need to be careful about
# what classes we use in defining
# objective functions in TMB.
#
# Success might indicate that we
# don't need to be careful, but
# there is always the possibility
# that we didn't choose an offending
# class.
# ------------------------------------------

# simulated incidence time series using the following parameters:
#   beta = 0.05
#   gamma = 0.01
#   nbdisp = 20
# and the following constants:
#   N = 1e7 -- total population size
#   I0 = exp(-3) ~ 0.04978707 -- initial number of infected individuals.
#                             -- now that I look at this, it doesn't
#                             -- look right, but I don't think this is an
#                             -- issue for the current goals
simulated_incidence =
c(21142, 33058, 17897, 31985, 26137, 25229, 29885, 37729, 39970,
  34085, 39473, 30872, 35766, 33358, 50120, 53363, 41875, 50854,
  31798, 47244, 46007, 39505, 61416, 62898, 57462, 59592, 43489,
  66236, 72622, 65582, 72186, 73235, 97367, 73209, 85605, 76294,
  110397, 111123, 124301, 67211, 117142, 80442, 76282, 119964,
  114966, 92815, 110600, 103442, 129273, 115845, 115017, 144605,
  93493, 115627, 89845, 124580, 134464, 128302, 117909, 142681,
  196050, 136264, 139951, 128954, 165404, 215344, 178069, 107529,
  244738, 170347, 174966, 234961, 197407, 132084, 148747, 125393,
  147786, 123326, 207731, 110553, 152299, 188208, 182393, 144416,
  151706, 92410, 147819, 81283, 68689, 96493, 92547, 142112, 121699,
  114169, 94884, 54743, 52095, 52281, 67706, 54149)

plot(simulated_incidence)

obj = MakeADFun(
  data = list(
    obs_incidence = simulated_incidence,
    N = 1e7,
    ## ----------------------------------------------------
    ## THIS IS A KEY BIT:
    ## ----------------------------------------------------
    ## - pass a sparse matrix in for beta+1, which is a scalar
    ## - this matrix will be modified on the C++ side,
    ##   and used in the computation of the likelihood
    ## ----------------------------------------------------
    beta_p_1 = as(Matrix::sparseMatrix(i = 1, j = 1, x = 0), "dgTMatrix")
  ),
  # start away from the true parameters, so that the
  # optimizer will take us back if it can
  parameters = list(
    log_beta = 0,
    log_gamma = 0,
    log_nbdisp = 0
  ),
  DLL = "SIR_try_to_break_tmb"
)

# optimize
obj$hessian <- TRUE
opt <- do.call("optim", obj)

# we converge, opt$convergence = 0
opt$convergence

# we get close to the true values
opt$par
log(c(beta = 0.05, gamma = 0.01, nbdisp = 20)) # true parameters

# not quite sure how to interpret these, but the
# gradient looks smallish -- can't remember how
# to scale gradients to assess convergence.
# pretty high correlation between log_beta and
# log_gamma, but this shows up whether we use
# log_beta directly in the likelihood calculations
# or via beta_p_1.
obj$fn(opt$par)
obj$gr(opt$par)
cov2cor(solve(obj$he(opt$par)))

# plot the fitted model
sims = obj$simulate()
lines(exp(sims$log_incidence), col = 'red', lwd = 3)

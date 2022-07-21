##########################################################
# Example of passing several vectors from R to C++ 
#     The number of vectors are flexible
#     These vectors don't have the same length
# Multiple matrices can be passed in the same way
##########################################################
library(TMB)
set.seed(1L)
n = 10

a = list(
  data = list(
    # These vectors (w or w/o names) are packed into a list
    Yx = list(Y = rnorm(5), x = rnorm(10), rnorm(3))
    #MatrixList = list(matrix(0, 1, 1), matrix(1, 2, 2)),
  ),
  param = list(
    a = as.numeric(0),
    b = as.numeric(0),
    logSigma = as.numeric(0)
  )
)
compile('strucs.cpp')
dyn.load(dynlib("strucs"))

a$data
print("===== before MakeADFun =====")
f = MakeADFun(data = a$data, parameters = a$param, DLL = 'strucs')
print("===== after MakeADFun =====")

#optim(f$par, f$fn)

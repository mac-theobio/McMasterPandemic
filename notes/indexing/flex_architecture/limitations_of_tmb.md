# Limitations of TMB

We have been unable to implement the following features in TMB.
* Calculate the dominant eigenvector of the Jacobian
* Round the state vector

These limitations appear to be related to unimplemented type casting in the automatic differentiation package that TMB leverages, CppAD.  Quote from Weiguang:
> The reason behind the failure is that TMB considers the cast as cast from type 'CppAD::AD<CppAD::AD<double> >' to type 'int' . CppAD::AD is a special type (different from typical types) in the sense that it not only stores a value of a variable but also keeps track of "symbol" of the variable through assignment and/or involvement in formula in order to do differentiation correctly. (edited) 


The RcppEigenAD package has an interesting example that does automatic differentiation of an eigenvector using CppAD, suggesting that it can be done.  But the plain example doesn't compile for me out of their box.

To avoid a potential similar issue I am avoiding the ceiling function in `make_delay_kernel` so that `max_len` just gets passed in through the `DATA` macro.

## notes on reimplementation

We want the model to be faster and more numerically stable. At present the log-likelihood surface seems very "rough", which we don't entirely understand but presume to be at least in part due to numerical instability. As a result, we are using *differential evolution* optimization to find the parameters; this means we have to evaluate the likelihood of thousands of parameter sets for each calibration, which means it takes 15-20 minutes or so.  (We can speed this up by running some of the evaluations in parallel, but we often want to calibrate lots of places at once, or try calibrating lots of variations of the same model, so this doesn't help too much.) What are our options for speeding things up?

* "C snippets": critical pieces of the model written in C, see [Morgan Kain's COVID interventions model](https://github.com/morgankain/COVID_interventions) or [CEID model](https://github.com/CEIDatUGA/COVID-stochastic-fitting). 
   * advantages: relatively straightforward
   * disadvanges: hard to debug
* `Rcpp`: write in C++
   * "sugar": lots of conveniences like vectorized operations, utility functions from R, matrix operations, etc..
   * huge documentation base about Rcpp (don't know if there are any COVID models written in Rcpp)
* `TMB`: template model builder, see [here](https://github.com/kaskr/adcomp) and [here](https://kaskr.github.io/adcomp/_book/Introduction.html). See [here](funcresp.cpp) for an "intermediate complexity" model under development from another project. I started to rough out an example [here](sim_funs.cpp)
   * calculates gradient and objective function simultaneously
   * uses the Laplace approximation to allow for integration across *random effects*, e.g. to account for process error or to fit several regions at once 
   * can be used with the `tmbstan` package to do Hamiltonian Monte Carlo
   
* need to start by implementing the equivalents of  `do_step`, `make_ratematrix`, `update_foi`, `make_state`, `run_sim` (we should eliminate `run_sim_range` as an intermediate step, instead constructing a fully filled-in timevars data frame of the changing parameters on every day of the simulation and passing this as data); we also need to compute the likelihood by comparing the results of the sim to the data

See also https://mrc-ide.github.io/squire/

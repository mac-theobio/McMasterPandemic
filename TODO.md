## to do

### short term

* epiestim: delay-convolve beta curve
* check off-by-one-day in breakpoints?

* `condense` option for including/not including ICU in H output
* add global options for condense: diff-deaths, hospital incl ICU, etc. ?

* fit to hospitalization only? or just leave out ICU?
* add testing intensity?
* MCMC
* arrange caching for MLi data (want to be able to replicate if data gets weird)
* make `get_r` work for r<1 ...
* distributions in `calibrate`
     * implement Poisson?
     * allow var-specific nb params
* `calibrate`/`mle_fun`: better error handling if dates/vars don't match
* `calibrate` substitute mle2 for optim?
* importance weights for `forecast_ensemble`?
* include hessian=TRUE by default
* process error version

* transition TODO list to GH issues??? (items marked Z have been posted)
* fix make rules
* Z write_params
* incorporate het/behavioural change: `beta*(S/N)^het_alpha*I`
* multi-run plots (aggregate + bind_rows)	
* improve documentation!

### minor/cosmetic

* add hospital admissions 'stock'
* Z translate label names for graphs ("D" -> "Dead", etc.)
* Z include start dates as metadata for params, states, runs ... 
* check `do_hazard` option and switch to this as default
* un-FRY `update_foi` function
* Z `write_params`, `read_params`; keep labels, unevaluated values, etc. as attributes, get start date as attribute. (Store as separate column, or commented row, or ... ???) [store params as JSON???]
* S3 method for params
    * Z pretty-printing (print rates as reciprocals?)
	
### longer term

* full stochastic version (now does obs stoch)
    * demographic: reulermultinom etc.
    * random-walk of parameters
* Bayesian/importance-sampling solution
* Erlang-izing transition matrix
* age structure
* spatial structure
* incorporate testing (see `testing_flow.md`)
* porting to other platforms for speed/latent variable handling 
     * TMB-ize ?
	 * generate enums of state names and params
	 * changing transition matrix in place will be *easier*

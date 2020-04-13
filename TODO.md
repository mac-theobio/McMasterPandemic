## to do

### short term

* implement Poisson in calibrate?
* allow var-specific nb params
* better error handling if dates/vars don't match
* roll in mle2???
* include hessian=TRUE by default

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

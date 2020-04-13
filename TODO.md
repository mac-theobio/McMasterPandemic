## to do

### short term

* encapsulate ensemble-forecast machinery
* implement Poisson in get_breaks_gen ?
* wrap new calibration machinery back into package
* better handling of multiple vars (auto-recognition of vars)
* better error handling if dates/vars don't match
* report non-convergence
* roll in mle2???

* transition TODO list to GH issues??? (items marked Z have been posted)
* fix make rules
* Z write_params
* incorporate het/behavioural change: `beta*(S/N)^het_alpha*I`
* fix time-varying machinery: check for off-by-one on dates?? (BMB)
* testing calibration (BMB)
* multi-parameter run interface (fix fitfun/predict/etc.)
    * make a predict method for parameter sets (assuming we can use E0 and start date only as part of params)
* multi-run plots (aggregate + bind_rows)	
* improve documentation!
* fitting procedure with time-varying parameters

### minor/cosmetic

* add hospital admissions 'stock'
* Z translate label names for graphs ("D" -> "Dead", etc.)
* Z include start dates as metadata for params, states, runs ... 
* check `do_hazard` option and switch to this as default
* un-FRY `update_foi` function
* Z `write_params`, `read_params`; keep labels, unevaluated values, etc. as attributes, get start date as attribute. (Store as separate column, or commented row, or ... ???) [store params as JSON???]
* S3 method for params
    * Z pretty-printing (print rates as reciprocals?)
	* Z transform/update method?

	
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

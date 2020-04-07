## to do

### short term

* transition TODO list to GH issues???
* fix make rules
* fix 
* write_params
* incorporate het/behavioural change: `beta*(S/N)^het_alpha*I`
* fix time-varying machinery: check for off-by-one on dates?? (BMB)
* testing calibration (BMB)
* multi-parameter run interface (fix fitfun/predict/etc.)
    * make a predict method for parameter sets (assuming we can use E0 and start date only as part of params)
* multi-run plots (aggregate + bind_rows)	
* get public H/D/ICU data for examples
* improve documentation!
* fitting procedure with time-varying parameters

### minor/cosmetic

* add hospital admissions 'stock'
* translate label names for graphs ("D" -> "Dead", etc.)
* include start dates as metadata for params, states, runs ... 
* check `do_hazard` option and switch to this as default
* carry labels as metadata with params
* un-FRY `update_foi` function
* `write_params`, `read_params`; keep labels, unevaluated values, etc. as attributes, get start date as attribute. (Store as separate column, or commented row, or ... ???) [store params as JSON???]
* expose 'safe' code in public-facing repo?
* S3 method for params
    * pretty-printing (print rates as reciprocals?)
	* transform/update method?


### longer term

* effective heterogeneity: alpha  = (S(t)/S(0))^alpha
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

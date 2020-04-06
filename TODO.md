## to do

### short term

* fix time-varying machinery: check for off-by-one on dates?? (BMB)
* testing calibration (BMB)

* allow internal, smaller time steps [TEST]
* test/fix calibration procedure
* multi-parameter run interface (fix fitfun/predict/etc.)
    * make a predict method for parameter sets (assuming we can use E0 and start date only as part of params)
* get public H/D/ICU data for examples
* improve documentation!
* fitting procedure with time-varying parameters

### minor/cosmetic

* add hospital admissions 'stock'
* translate label names for graphs ("D" -> "Dead", etc.)
* include start dates as attributes of params
* include parameters and other metadata as an attribute in returned object
* check `do_hazard` option and switch to this as default
* enable dt<1 for testing
* un-FRY `update_foi` function
* `write_params`, `read_params`; keep labels, unevaluated values, etc. as attributes, get start date as attribute. (Store as separate column, or commented row, or ... ???) [store params as JSON???]
* expose 'safe' code in public-facing repo?
* S3 method for params
    * pretty-printing (print rates as reciprocals?)
	* transform/update method?


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

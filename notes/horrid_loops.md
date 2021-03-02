An attempt to document what is currently going wrong with `sim_funs. Interaction between `make_state()`, testify machinery, eigenvector machinery for generating initial conditions ... We switched two defaults: (1) start from eigenvector, (2) use testify machinery if the current parameters have a testing-rate parameter >0.

I think the fundamental design problems are that 

* we have a tall stack of functions (`calibrate_comb()` → `calibrate()` → `run_sim_loglin()` → `run_sim()` → `do_step()`) and we are trying to hard to let users call any of those functions and have them just "do what you mean", i.e. flexible/clever defaults
* there is a bit of recursion going on: `make_ratemat()` → `make_state()` (default: may automatically try to compute eigvecs and/or testify) → `make_ratemat()` (if eigvecs desired, to set up the exponential simulation)

What happens in `make_state(params=pp)`?

* figure out `testify` (`params` specified and has testing>0?)
* set base state names (no testing/age etc.)
* `expand_stateval_testing(state, method="untested")`
   * call `get_evec(params, testify=testify)`
      * call `rExp()`
	     * call `make_state(N=1,E0=1e-5,use_eigvec=FALSE, testify=FALSE)`
		 * call `make_ratemat(state=state)` (no testify)
		 * testify `M`
		 * `expand_stateval_testing()`
		 * `run_sim_range(params, state, M)`
   * assign `evec` values back

* sort out starting/ending dates etc. (harmless)
* if `state` not specified, call `make_state` **with testify=FALSE**
    * call `get_evec()`
	   * call `rExp()`
	       * call `make_state(N=1,E0=1e-5, use_eigvec=FALSE)`
		   * call `make_ratemat(
		          
* call `make_ratemat`

    * this uses `state`, so will call `make_state(..., params=params)` [except for infection terms, which are determined by `update_foi`, `make_state()` only needs to know about the dimensions and names of the state vector: could simplify?]
	

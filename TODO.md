## to do

### short term

* allow passing NLLvec to `pop_pred_samp`
* per-variable observation noise
* ?? fix "predict three times"
* test spline variations: number of spline parameters? effects of var choice on knot placement? effects of (strong) penalization on spline coeffs?

* document !!!!
* document Rt calculation (add to McMasterReport.pdf) etc.
* allow choice of Wald vs DEoptim cov matrices in forecast_ensemble; allow control of `pop_pred_samp`
* document invlink_trans better
* secret plots; writeup for secret

## important

* save fluxes
* calibrate uncertainty
    * DEoptim pop/imp wt experiments
	* full Bayes?
    * pop_pred_samp
* match up with EpiEstim
* LTCF: data, model compartments
* TMB core
* testing flow

## debugging

* should log_mu be logit_mu instead?
* debug update/parallel DE_cores

### strategic/architectural

* allow browser()/debugging in mle_fun (circumvent mle2 problems)?
* condense shouldn't be so hard-coded; refactor. Enforce conventions for var names (I vs ICU)?
* think about condensation/diff/cum (cumRep, hosp admission.  Can cumulate during or post-run; need something that makes obs error correct and is efficient)
* Erlangization/chain trick
* TMB/Stan implementations
    * compartmental
    * renewal
* testing (distributions, compartmental)	
	
### technical dept/consistency

* change/document that parameters should *not* contain underscores (or use a different delimiter for invlink prefixes, e.g. .. or __ or |
* rename fix_pars to fix_params?
* print start, end dates in print.fit_pansim method
* arrange consistent factor ordering/colour palette throughout
* simulations with per-variable obs error, ¿process error?
* hospital admissions
* pipeline for copying up-to-date calibrations and data into package
* priors on r/R0/etc.?
* MCMC
* carry along last date as metadata in calibrations (what data were used?)
* process error
   * implemented but not well-documented
   * overdispersion (beyond demog stoch) is global, not per-transition or per-state variable ... 
   * allowing per-trans overdisp means thinking about including vectors in a parameter list - what if anything will this break?
* R(t)
   * delay-convolve beta curve
   * consider other estimation machinery?
* check hospital-only calibration (don't try to calibrate)
* Makiness: 
    * depend on package version
	* ont_all from clean, not calib?
	* dotdir

* add testing intensity?
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

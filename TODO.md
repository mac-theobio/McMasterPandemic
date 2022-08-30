## ROADMAP

### Prioritization

![[Roadmap prioritization diagram on GitHub.](https://github.com/mac-theobio/McMasterPandemic/blob/master/assets/roadmap.svg)](assets/roadmap.svg)

### Descriptions of the Items

* Arithmetic Engine
  * Main existing pain-points:
    * Interface for inputting rate expressions is not intuitive in many (but not all) cases
    * Not all operations that a user would want to include in a rate expression are possible (e.g. polynomials with more than one term)
  * Details:
    * Ensure that operations can happen in any order (e.g. currently powers must be taken after sums of state variables)
    * User can just use the formula interface without needing to consider storing intermediate values
    * Use 'condensation operations' like convolution and lagged differencing more generally
    * Templatize the C++ code so that the user can compile it themselves with custom operations (e.g. "for some reason I just need the digamma function")
* Numerical Robustness
  * Main existing pain-points:
    * Users get no help on how to diagnose lack of convergence
    * Users get warnings on negative state variables and rates, but it is not clear how to figure out why and how to stop it
  * Details:
    * Currently `convergence_info`, `opt_obj`, and `opt_par` provide basic access to the underlying optimization objects, which should be enough for people to construct their own convergence diagnostics
    * The first step to making this easier is to finish [this](https://canmod.github.io/macpan-book/convergence.html)
    * Then we might consider wrapping this stuff up
    * Provide warnings for common scenarios (e.g. I've noticed that the default `bbmle` settings can get stuck on parameters that cause flat dynamics for simple models -- not sure how important this is for more complex models)
    * Recommended and wise default priors
    * Restricting state variables to be positive is difficult given that our formulation is not on the log scale
    * We need better documentation on how to run simulations to diagnose negativity
* Inference on Parameters
  * Main existing pain-points:
    * Not clear to users how to get confidence intervals for parameters
  * Details:
    * Can get confidence intervals from `bbmle::confint.mle2`, `TMB::sdreport`, `rstan::summary`
    * Should at least provide examples in the guide, and then wrap them if it makes sense
* Log-Linear Models
  * Main existing pain-points:
    * The current piece-wise model of time-variation limits modelling flexibility
  * Details:
    * More general parameter time-variation
    * Any parameter is a link-transformed linear matrix equation
    * We should add terms for dependence on state variables
    * How to extend semi-parametric log-linear models into forecasting periods?
* Process Error (Sim ProcErr and Calib ProcErr)
  * Main existing pain-points:
    * Process error has important dynamical consequences for many models, but is not currently possible
  * Details:
    * Should be easy to add to simulations, but it will be a little more challenging to include with calibrations
    * This will need to plug in nicely with log-linear models, where we can decompose the model matrix into a fixed-effect and random-effect component
    * The random-effect component will naturally make a fairly wide range of process error models available including 
* Model Products
  * Main existing pain-points:
    * Model definition might be currently harder than it needs to be, given the lack of a formal but straightforward interface for combining sub-models 
  * Details:
    * Model Specification as Products of Sub-Models
    * Multiply sub-models (e.g. `SEIR * vax_status`)
    * Currently possible in a lightweight way using string manipulation of names, so it is possible that we just need a better-documented set of examples
    * Variation of parameters over sub-models
      * i.e. like params_timevar but over other model categories 
      * e.g. transmission rate varies over vaccination status
    * We have [detailed notes that include a proposal](https://hackmd.io/@stevencarlislewalker/B1hxAvxBc)
* ODE Solvers
  * Main existing pain-points:
    * Within-time-step variation can have important dynamical consequences in some cases, but this is not currently possible
  * Details:
    * Runge-Kutta-4 should be straightforward
    * Is the boxcar model related to this?
* R0, Gbar, r
  * Main existing pain-points:
    * Applying default parameters to calibration/forecasting scenarios can lead to bad results, and fix-par functionality provides some protection
  * Details:
    * R0, Gbar, r for general models (next-generation matrix?)
    * Priors on these aggregate summarizing parameters
* Generic Time Step Units
  * Main existing pain-points:
    * The natural time-step is not always one day, but this is currently a hard requirement
  * Details:
    * Each step always means one day
    * But more often we only have data at less granular time scales
* Flows without per-capita Rates
  * Main existing pain-points:
    * Some flows that we encounter are unnatural or impossible to specify as per capita rates (e.g. Vital dynamics with waning;  vaccination administration)
  * Details:
    * e.g. vaccination
    * e.g. for convenience when the per-captia denominator is the sum of multiple boxes as it often is with birth rates
      * This is a convenience issue for birth rates but becomes a necessity when both birth and waning are in the model, because then waning needs to have an associated outflow but not birth -- and the engine has no ability to restrict outflow from some states back to S but not from others
* MCMC
  * Main existing pain-points:
    * I'm not sure that there are any existing pain-points
  * Details:
    * tmbstan vs rstan engine as an alternative to tmb
    * Probably should start with tmbstan as we already have an prototype
* Uneven Lagged Diff
  * Main existing pain-points:
    * Reports are not always provided every `n` days, and so the lagged differencing needs to account for variable `n`
* Technical Debt
  * Modularization & Software Stability
    * Main existing pain-points:
      * Too many parts of the codebase need to be modified when adding new features
      * The Clean up Existing Features aspect of technical debt grows faster when Modularization is insufficient
    * Details:
      * No new 'user' functionality
      * After all I've learned, redesign the codebase so that it is easier to extend in the future
      * Big value for adoption and community growth
      * C++ modularity versus R-side modularity -- which is more important?
      * Design principles
        * Single responsibility -- each function should have one reason to change
        * Open-closed -- functions should be open for extension but closed for modification
        * Substitution -- alternative versions of a method should have the same input and output types, but the computation differs
        * Interface segregation -- interfaces should not depend on methods that they do not use, and therefore we should tend to have many use-case focused interfaces rather than few general-purpose interfaces
        * Dependency inversion -- interfaces should depend on abstract function definitions, not concrete implementations
  * Clean up Existing Features
    * Main existing pain-points:
      * The list is growing, which increases developer stress
    * Details:
      * Enforce Ordering of Model Definition Functions
      * Ensemble Forecasting
        * Choose any combination of parameter uncertainty, observation error, and process error (when process error is possible)
        * Flexibly omit error for particular time-series or omit uncertainty for certain parameters
        * Ensemble forecasting should run faster in loops within simulation macros on the C++ side, as opposed to generating samples from the distribution of parameters on the R side and then passing back to C++
        * Smoothing out the empirical quantiles so that forecast envelopes are not as bumpy
      * Symbolic Matrix Algebra Engine
        * Poorly documented
        * If we do a cleanup of the Arithmetic Engine then the matrix stuff should be done at that time, it we do not then the matrix stuff needs to be done sooner
      * Make State Cleanup
        * Option to not use eigenvector and just populate the initial state with elements of the parameter vector
        * More intuitive interface, or at least detailed documentation
      * Condensation
        * Gaussian convolution -- make sure that it is robust when sd is large (set a threshold for an error -- in general we need this for all convolutions)
          * Mass loss problem: the q-vector that we choose will cause us to lose mass at the tails
            * Could renormalize (probably should)
            * Take the difference of the cumulative convolution distribution of the two end-points, and if this is above a threshold then throw an error/warn/fix??
            * We could also use this approach to check for negativity
        * Cumulative sums in condensation on the C++ side
      * Ability to fit prior hyperparameters
      * Outflow
        * The auto-outflow approach is confusing, because people rightly expect add_outflow to _not_ mean subtract outflow but it can mean exactly that depending on how the model is set up
  * Guides
    * Main existing pain-points:
      * We have a [draft but reasonably helpful guide for users](https://canmod.github.io/macpan-book/), but the [contributors](https://github.com/mac-theobio/McMasterPandemic/blob/master/CONTRIBUTING.md) guide is totally insufficient
    * Details:
      * Part of the problem is that the previous technical debt items demotivate me from finishing the above



### Old TODO List -- Still Useful!

#### very short term

Naming conventions for value vector?
troubleshoot test_calibrate
more docs for new params_timevar structure
issue/document testLevel stuff in README

update method for calibrations

names?  Symbol.Date
values -> time_params

#### short term

* **REAL CONFIDENCE INTERVALS!!**
   * we need to get importance weighting going etc. (or MCMC, but imp weighting should work if anything will)

* instabilities etc.
     * rewrite as gradient
     * Runge-Kutta steps?
	 * fancy flows/transformations
	 * why are hospitalizations bouncy?
	 
* simulate from simple model
* simulate from spline model

* test spline variations: number of spline parameters? effects of var choice on knot placement? effects of (strong) penalization on spline coeffs?

* document !!!!
* document Rt calculation (add to McMasterReport.pdf) etc.
* allow choice of Wald vs DEoptim cov matrices in forecast_ensemble; allow control of `pop_pred_samp`
* document invlink_trans better
* secret plots; writeup for secret

### important

* save fluxes
* calibrate uncertainty
    * DEoptim pop/imp wt experiments
	* full Bayes?
    * pop_pred_samp
* match up with EpiEstim
* LTCF: data, model compartments
* TMB core
* testing flow

### debugging

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

# McMasterPandemic Roadmap

* Arithmetic Engine
  * Main existing pain-points:
    * Interface for inputting rate expressions is not intuitive in many (but not all) cases
    * Not all operations that a user would want to include in a rate expression are possible (e.g. polynomials with more than one term)
  * Details:
    * Ensure that operations can happen in any order (e.g. currently powers must be taken after sums of state variables)
    * User can just use the formula interface without needing to consider storing intermediate values
    * Use 'condensation operations' like convolution and lagged differencing more generally
    * Templatize the C++ code so that the user can compile it themselves with custom operations (e.g. "for some reason I just need the digamma function")
* Convergence Diagnostics and Numerical Robustness
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
* Joint Observation Errors
  * Main existing pain-points:
    * Variation due to observation error does not respect the constraint that fluctuations should average out over the compartments
  * Details:
    * This seems really hard in the face of reporting delays, and may is so hard that we can ignore it?
* Process Error
  * Main existing pain-points:
    * Process error has important dynamical consequences for many models, but is not currently possible
  * Details:
    * Should be easy to add to simulations, but it will be a little more challenging to include with calibrations
    * This will need to plug in nicely with log-linear models, where we can decompose the model matrix into a fixed-effect and random-effect component
    * The random-effect component will naturally make a fairly wide range of process error models available including 
* Model Specification as Products of Sub-Models
  * Main existing pain-points:
    * Model definition might be currently harder than it needs to be, given the lack of a formal but straightforward interface for combining sub-models 
  * Details:
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
* Fix Pars
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
* Uneven Lagged Differencing
  * Main existing pain-points:
    * Reports are not always provided every `n` days, and so the lagged differencing needs to account for variable `n`
* Technical Debt
  * Modularization
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
      * Symbolic Matrix Algebra Engine
        * Poorly documented
        * If we do a cleanup of the Arithmetic Engine then the matrix stuff should be done at that time, it we do not then the matrix stuff needs to be done sooner
      * Make State Cleanup
        * Option to not use eigenvector and just populate the initial state with elements of the parameter vector
        * More intuitive interface, or at least detailed documentation
      * Ensemble Forecasting
        * Choose any combination of parameter uncertainty, observation error, and process error (when process error is possible)
        * Flexibly omit error for particular time-series or omit uncertainty for certain parameters
        * Ensemble forecasting should run faster in loops within simulation macros on the C++ side, as opposed to generating samples from the distribution of parameters on the R side and then passing back to C++
        * Smoothing out the empirical quantiles so that forecast envelopes are not as bumpy
      * Condensation
        * Cumulative sums in condensation on the C++ side
        * Gaussian convolution -- make sure that it is robust when sd is large (set a threshold for an error -- in general we need this for all convolutions)
        * Mass loss problem: the q-vector that we choose will cause us to lose mass at the tails
          * Could renormalize (probably should)
          * Take the difference of the cumulative convolution distribution of the two end-points, and if this is above a threshold then throw an error/warn/fix??
          * We could also use this approach to check for negativity
      * Ability to fit prior hyperparameters
  * Guide for Contributors
    * Main existing pain-points:
      * We have a [draft but reasonably helpful guide for users](https://canmod.github.io/macpan-book/), but the contributors guide is totally insufficient
    * Details:
      * Part of the problem is that the previous technical debt items demotivate me from finishing the above

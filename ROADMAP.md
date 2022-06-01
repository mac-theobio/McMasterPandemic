# McMasterPandemic Roadmap

* Arithmetic Engine
  * Ensure that operations can happen in any order (e.g. currently powers must be taken after sums of state variables)
  * User can just use the formula interface without needing to consider storing intermediate values
  * Use 'condensation operations' like convolution and lagged differencing more generally
  * Templatize the C++ code so that the user can compile it themselves with custom operations (e.g. "for some reason I just need the digamma function")
* Log-Linear Models
  * More general parameter time-variation
  * How to extend semi-parametric log-linear models into forecasting periods?
* Joint Observation Errors
  * Variation needs to respect that the total population size is not affected by variation -- this seems really hard in the face of reporting delays, and may is so hard that we can ignore it?
* Process Error
  * Initial step would be to treat params_timevar values as random effects
  * This will need to plug in nicely with log-linear models
* Model Specification as Products of Sub-Models
  * Multiply sub-models (e.g. `SEIR * vax_status`)
  * Currently possible in a lightweight way using string manipulation of names, so it is possible that we just need a better-documented set of examples
  * Variation of parameters over sub-models
    - i.e. like params_timevar but over other model categories 
    - e.g. transmission rate varies over vaccination status
* Runge-Kutta-4 ODE Solvers
* Fix Pars
  * R0, Gbar, r for general models (next-generation matrix?)
  * Priors on these aggregate summarizing parameters
* Generic Time Step Units
  * Each step always means one day
  * But more often we only have data at less granular time scales
* Flows without per-capita Rates
  * e.g. vaccination
  * e.g. for convenience when the per-captia denominator is the sum of multiple boxes as it often is with birth rates
    * This is a convenience issue for birth rates but becomes a necessity when both birth and waning are in the model, because then waning needs to have an associated outflow but not birth -- and the engine has no ability to restrict outflow from some states back to S but not from others
* MCMC
  * tmbstan vs rstan engine as an alternative to tmb
  * Probably should start with tmbstan as we already have an prototype
* Uneven Lagged Differencing
* Restricting State Variable Values (e.g. positivity)
  * Need to respect differentiability
* Technical Debt
  * Modularization
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
    * Symbolic Matrix Algebra Engine
    * Make State Cleanup
      * Option to not use eigenvector
    * More Flexible Ensembles
      * Choose any combination of parameter uncertainty, observation error, and process error
      * Flexibly omit error for particular time-series or omit uncertainty for certain parameters
    * Condensation
      * Cumulative sums in condensation on the C++ side
      * Gaussian convolution -- make sure that it is robust when sd is large (set a threshold for an error -- in general we need this for all convolutions)
      * Mass loss problem: the q-vector that we choose will cause us to lose mass at the tails
        * Could renormalize (probably should)
        * Take the difference of the cumulative convolution distribution of the two end-points, and if this is above a threshold then throw an error/warn/fix??
        * We could also use this approach to check for negativity
    * Ability to fit prior hyperparameters

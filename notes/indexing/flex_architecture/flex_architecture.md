## Concepts

* A `compartmental model` has the following sub-models
  * State structure
  * Default parameters (see formal definition of parameters below)
  * Intermediate variables (e.g. sum of infective compartments)
  * Symbolic vectors and matrices (e.g. vaccine transmission reduction matrix)
  * State initialization steps (e.g. eigen-vector of a state-reduced jacobian -- I know ... this is unclear)
  * Rates (e.g. per-capita rates of state transition)
  * Outflow (what states transitions remove individuals from the origin state?)
  * Post processing steps (e.g. convolution of incidence into reports)
* A `simulation model` = `compartmental model` + `simulation engine` + `simulation schedule`
* Distinction between parameters and settings
  * Parameters are numerical variables (floats or doubles, but not integers) whose values can be changed after the `simulation model` is rendered
  * Settings are variables that establish the structure of the model, and whose values cannot be changed after the `simulation_model` is rendered
* There are different types of parameters for different sub-models
  * Compartmental Model Parameters
    * State initialization (e.g. initial number of exposed individuals)
    * Simulation step (e.g. transmission rate)
    * Post processing (e.g. convolution delay mean)
  * Simulation Schedule Parameters
    * Time variation (e.g. breakpoint multiplier)
    * Observational stochasticity (e.g. negative binomial dispersion parameter)
    * Observational data
    * Error distribution (e.g. negative binomial dispersion parameter)
  * Epidemiological statistics (e.g. R0)


## Model Structure Information



## Sub-Model Short Names and Associated Abbreviated Prefixes

* Compartmental -- `box_` -- maybe `bx_`
  * Symbolic -- `sym_`
  * State -- `state_` -- maybe `st_`
  * Simulation -- `sim_`
  * Post -- `post_` -- maybe `pp_` for 'post processing'
* Temporal -- `temp_` -- maybe `tl_`
  * Bounds -- `bounds_` -- maybe `rng_` for 'time range'
  * Special -- `spec_` -- maybe `sd_` for 'special day'
  * Exogenous -- `exog_` -- maybe `ex_`
  * Observed -- `obs_`


## Examples of Parameters and Settings in Different Sub-Models

* Compartmental Model
  * State initialization sub-model
    * Settings
      * Character vector of state names
      * Values for initial state vector (optional if the initial state is the same for all parameter sets)
      * Settings for how to use the eigenvector of the initial Jacobian to populate the initial state vector (e.g. what states are exposed?)
    * Parameters
      * Initial number of exposed individuals
      * A set of values that is to be copied over to the state vector to initialize it
  * Simulation step sub-model
    * Settings
    * Parameters
      * Transmission rate
      * Recovery rate
      * Population size
      * Layer-specific rates (e.g. unvaxed vs one dose vs two doses etc ...)
  * Post processing sub-model
    * Settings
    * Parameters
      * Convolution delay kernel mean
      * Negative binomial dispersion parameter for observational stochasticity
* Temporal Localization Model
  * Simulation bounds sub-model
    * Settings
      * Start/end date
      * Start date offset for calibrations
      * End date offset for forecasts
    * Parameters
  * Special day sub-model
    * Settings
      * Declare weekend effects
    * Parameters
  * Exogenous time variation sub-model
    * Settings
    * Parameters
      * Breakpoint multiplier
      * Smoothing spline coefficient
  * Data comparison sub-model
    * Settings
      * Observed dataset
    * Parameters
      * Dispersion parameter of the error distribution



## Data Structure Model

Hypothesis: the model definitions can be thought of as a series of tables that are related by various unique keys, much like a database.

* Compartmental Model
  * States
    * state -- primary key
    * structural layers or dimensions (e.g. vaccination category)
  * Parameters
    * parameter -- primary key
    * structural layers or dimensions (e.g. vaccination category)
    * default values
  * Summaries
    * name -- primary key
    * type -- sum, intermediate, lag-n difference, convolution
    * expression -- expression defining the summary in terms of states, parameters, and other summaries
  * Rates
    * from -- foreign key to match with State; primary key together with To
    * to -- foreign key to match with State; primary key together with From
    * expression -- expression defining the rate in terms of states, parameters, and/or summaries
    * outflow -- true/false whether these rates induce outflows as well as inflows
* 

Do we split up different types of parameters?

## Sub-Model Compatibility

* Symbolic-State: nothing
* Symbolic-Simulation: some names in formulae must match matrix variable and dimension names
* Symbolic-Post: some names in 
* State-Simulation:
* State-Post:
* Bounds-Special:
* Bounds-Exogenous:
* Bounds-Observed:
* Special-Exogenous:
* Special-Observed:
* Exogenous-Observed:
* Compartment-Temporal: 

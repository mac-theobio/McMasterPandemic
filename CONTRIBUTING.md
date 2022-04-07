# Contributing to McMasterPandemic

Thank you for contributing to our project.

## Architecture

### Current Architecture

The simplified current architecture is described by the following diagram.

![](assets/architecture.svg)

Each box is a layer of functions that call functions in the layers below them.

Model definition functions are used to create and manipulate parameter and state-variable objects that define model structure and defaults values.  These objects have S3 classes `params_pansim` and `state_pansim`.  Common functions in this layer are `read_params`, `make_state`, and `expand_state_vax`.

The forecasting and simulation functions allow one to ... `do_step`, `run_sim`, `forecast_sim`, and others that are typically not used by users (e.g. `run_sim_range`, ...)

### Target Architecture Proposal

![](assets/target_architecture.svg)

### Engine-Client Separation

* Currently one client -- the R client -- and hopefully/probably we will only ever have one client
* Currently two engines
  * Original R
  * TMB


* Stuff that goes in the engine:
  * Indices
  * Simulation steps
  * Simulation output summarization (e.g. convolutions, differences)
  * Objective function computation
* Stuff that goes in the client:
  * Names
  * Names of states
  * Names of parameters (e.g. transmission rate)
  * Names of loss functions (e.g. negative binomial, gaussian)
  * Topology of flow among states

#### NOTES

* What about summary methods for parameters?
  * Use case -- compute R0/r/Gbar for an _ensemble_ of mcmc simulations of parameters
  * In general we might need a many-to-one relationship between parameter sets and compartmental models
* In general we need to think about ensembles, sampling from parameter distributions, importance sampling, mcmc, etc...


## Engines



## Maintaining C++ Code

We try to be disciplined in how we maintain the C++ code


## Comparing R and TMB Engines

TODO: describe global options, `tmb_mode`, `R_mode`, etc ...

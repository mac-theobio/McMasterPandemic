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

We try to be disciplined in how we maintain the `C++` code, in that the correctness of this code can in principle be judged against [this specification document](https://canmod.net/misc/flex_specs).  These specifications are [versioned](https://canmod.net/misc/flex_specs#versioning-and-lifecycle) with the reference implementation for each version located in [inst/tmb](https://github.com/mac-theobio/McMasterPandemic/tree/tmb/inst/tmb). The currently recommended version is saved in [inst/tmb/recommended_spec_version](https://github.com/mac-theobio/McMasterPandemic/tree/tmb/inst/tmb/recommended_spec_version).

The package itself can only use one of these reference implementations at a time. The default spec version can be set by copying the appropriate versioned `cpp` file into `src/McMasterPandemic.cpp`, and this can be accomplished using the following [`make`](https://www.gnu.org/software/make/) rule.
```
make src/McMasterPandemic.cpp
```
This `make` rule copies the currently recommended version into `src/McMasterPandemic.cpp`.

One may also switch the spec version being used during an `R` session by using `set_spec_version` function. For example, if executed from the project root directory the following line will update the spec version and associated `C++` code being used.
```
set_spec_version('0.1.0', 'inst/tmb')
```
Note that this will also likely trigger compilation of the `C++` code unless the objects, shared objects, and/or DLLs can be found in the same directory as the associated `.cpp` file or some method of caching has been established.  One may find the spec version that they are using at any time with the following function.
```
spec_version()
```

On the `R` side one may also test for a particular spec version using the following family of functions: `spec_ver_eq`, `spec_ver_gt`, `spec_ver_lt`, `spec_ver_btwn`. Similarly, one may assert a certain spec version using the `spec_check` and `feature_check` functions.


## Comparing R and TMB Engines

TODO: describe global options, `tmb_mode`, `R_mode`, etc ...

TODO: describe `compare_*` family of functions ...on 


## The `flexmodel` Class

TODO: point to a help file (not yet written) that describes the `flexmodel` class

TODO: describe the `get_*` family of functions and how they should be used to create functions like `rate_summary`


## Global Options

To identify McMasterPandemic-specific global options on the `R` side we prefix their names with `MP_`. For example, the `MP_use_state_rounding` option is used to indicate whether the R engine should ever round the state variable. The defaults for all McMasterPandemic-specific options should be set in [R/zzz.R](https://github.com/mac-theobio/McMasterPandemic/blob/tmb/R/zzz.R).

## Package Tests



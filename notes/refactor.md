How do we go about modifying/refactoring the code to allow it to more flexibly/generally handle different models of time-varying parameters?

At the moment the calling sequence is: 

```
calibrate → mle_fun → forecast_sim → run_sim_breaks → run_sim → run_sim_range
```

(`predict` also calls `forecast_sim`, directly and via `forecast_ensemble`)

- `forecast_sim`: takes parameters
     - `p` (numeric parameter vector)
	 - `opt_pars` (list: starting values *and structure* for "re-listing")
	 - `start_date`, `end_date`
	 - `break_dates`
inverse-link transforms parameters, re-lists them into the format of `opt_pars`, then calls `run_sim_break`; then condenses/aggregates/pivots results to match data format 
	 
 - `run_sim_break`: thin wrapper for `run_sim`: converts `break_dates` plus relative beta information (`rel_beta0`) into a data frame of points at which parameters change, then passes that info plus parameters to `run_sim`
 
- `run_sim`: constructs rate matrix; processes information about when to turn on stochasticity. Loops over change points, calling `run_sim_range` repeatedly

- `run_sim_range`: assumes constant *per capita* rates (except for infection, which is constan
 
 

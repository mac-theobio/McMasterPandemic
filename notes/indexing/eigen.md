# non-vax case

* pass `params` -- are they modified?
* `do_hazard` is `FALSE` by default
* `state` is coming in `NULL`

* set N = 1
* call `make_state` with `N = 1` and `E0 = 1e-05`, `use_eigvec = FALSE`
  * Inside `make_state`
  * initialize state to the basic set of variables `"S"    "E"    "Ia"   "Ip"   "Im"   "Is"   "H"    "H2"   "ICUs" "ICUd" "D"    "R"` all set to zero
  * Set all `E` states to `E0`
  * Set all `S` states to `N - sum(E0)`, which I think really is just `N - E0`
  * Return `state`
* Make the `ratemat` with this `state`
* call `run_sim_range` with `do_exponential = TRUE` and `testwt_scale = "none"`
  * Inside `run_sim_range`
  * Initialize `foi` columns to `NA`
  * Initialize `res` matrix to `NA`
  * Loop over iterations
    * Call `do_step`
      * Inside `do_step`
      * `outflow[S_states] = rowSums(flows[S_pos, S_pos])`
      * `outflow[non_S_states] = rowSums(flows[notS_pos, p_states])`
* 
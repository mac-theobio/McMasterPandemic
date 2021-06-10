
`do_step_cpp` should take inputs:

- `state` (numeric vector)
- `ratemat` (**sparse** square matrix, dims match `length(state)`
- ? flags for `do_hazard`, `stoch_proc`, `do_exponential`, `testwt_scale` (to start, hardcode default behaviour TRUE, TRUE, FALSE, "N")
- indices of parallel accumulators?
- indices/information for estimating foi (definitely), testing & vacc (maybe??)

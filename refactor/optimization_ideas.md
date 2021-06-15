## optimization_ideas.md

We should put our ideas for speeding up the package/calibration here. Additionally, we should estimate the potential payoffs of implementing each idea.
* Use online calibration (estimated 50-75% speedup, from 3-4 hours instead of 7-14 hours). However, we still should do full calibrations once every 2 weeks or so. 
* run calibration on Graham or other supercomputer cluster (extreme speedup especially if given high priority)
* add_updated_vaxrate() copies the entire rate matrix, makes changes to it, then returns the updated rate matrix. I’m sure we can speed this up using some indexing tricks Ben mentioned in the last meeting. (up to 50% speedup, but most likely less than 25%)
  *   We'd have to think about it carefully, but one **Very** quick way to achieve pass-by-reference semantics in R (rather than the usual copy-on-modify) is to nest the item you want to update in place inside an environment.
* Reimplementing do_step in Rcpp (??? speedup)
* only use regex once, and save the resulting variable positions as an object attribute (estimated 15% speedup)
* replace regexes with startsWith (done, 3% speedup)
* optimize Kronecker product function in vaxify model (done, 5-10% speedup)
* replacing columnwise vector-matrix multiplication with M*v (done, 5-10% speedup)
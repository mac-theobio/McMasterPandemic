# Potentially Relevant Changes

* https://github.com/mac-theobio/McMasterPandemic/commit/b3973afe2ad859254adcea1d700d08d569c4e9f9
  * Expandify
  * make_vaccination_model
* https://github.com/mac-theobio/McMasterPandemic/commit/d7cbbed680896ae03d9e8ea89a73954846659211
  * Common multiples
    * get_var_indx, get_var_counts, get_var_modifiers
  * Eigenvector tolerance
  * Spec version bump to 0.1.2
  * (NO -- never used) `#include <cppad/local/cond_exp.hpp>`
* https://github.com/mac-theobio/McMasterPandemic/commit/36d5bcd5602b6a0b7cea1a5f5f7893a039f50061
  * `getOption("MP_tmb_models_match_r")`
* https://github.com/mac-theobio/McMasterPandemic/commit/dc2b2ccb1238c040e5a289d029d4235f66b1feb8
  * winstretch -- almost surely not the issue
* https://github.com/mac-theobio/McMasterPandemic/commit/b0f1c6c5082f35f07bbf50522a5a0cc9bc54a498
  * `nt <- round(gg[["Gbar"]] * winstretch)` -- this rounding _could_ be the issue, but probably not


# Hypotheses

* make_vaccination_model change
  * this is the most likely issue
  * can be tested by r_tmb_comparable()
* evec tolerance bug
  * this seems like the most likely to me
  * causes difference in initial state
  * might be able to prove a difference by comparing initial TMB-computed state on both branches
* factrs
  * seems very unlikely given that `make_vaccination_model` does not include the `factr` feature -- this has only been used for `threedosewane`
  * also this feature has been tested reasonably well -- although never in production
* expandify
  * likely not the issue, unless `expand_params_variant` is called by the pipeline
  * but even then ...
* winstretch rounding
  * perhaps this causes differences in fix_pars, but not sure
  * can mike test this?



single arrow from R_fully_vax back to S_fully_vax -- 
# Understanding `make_beta` Under Two-Dose Vaccination Model

* `Icats` are the names of the infected states -- `Ia, Ip, Im, Is`
* `Icat_prop_vec` is `Ca, Cp, (1 - iso_m) * Cm, (1 - iso_s) * Cs)`
* `beta_0` is `beta0 * 1/N * Icat_prop_vec`
* `vax_cat` is a list of state variable suffixes that identify vaccination status `unvax, vaxdose1, vaxprotect1, vaxdose2, vaxprotect2` for a two-dose model -- every state variable is categorized in this way
* `vax_trans_red` is a 5-by-5 matrix (corresponding to the `vax_cat` status')
  * set the `vax_cat == vaxprotect1` row to `1 - vax_efficacy_dose1`
  * set the `vax_cat == vaxdose2` row to `1 - vax_efficacy_dose1`
  * set the `vax_cat == vaxprotect2` row to `1 - vax_efficacy_dose2`
  * all other rows at `1` so far
* `beta_0` gets updated to `kron(vax_trans_red, t(beta_0))`
* `beta` is a matrix with 5 rows (same as `beta_0`) and `n_state` columns -- fill the corresponding `beta_0` values in `beta`, leaving the remaining columns all zero (so we can leave these out in a more efficient implementation)
* ****

* `testcats` are ??


# Understanding `add_updated_vaxrate` Under Two-Dose Model

* Compute the `vax_rate` list with one scalar element per dose
* Set the dimensions of the problem
  * `epi_states` -- `"S"    "E"    "Ia"   "Ip"   "Im"   "Is"   "H"    "H2"   "ICUs" "ICUd" "D" "R"    "X"    "V"`
  * `asymp_cat` -- `"S"  "E"  "Ia" "Ip" "R"`
  * `vax_cat` -- `"unvax" "vaxdose1" "vaxprotect1" "vaxdose2" "vaxprotect2"`
* Initialize two identical block matrices (`vax_block_dose1` and `vax_block_dose2`) with all zeros, with rows and columns each corresponding to one of the `epi_states`
* Fill in `vax_rate$dose_n` into each of the following elements of `vax_block_dose_n` (where `_n` is either 1 or 2)
  * Diagonal elements corresponding to each `asymp_cat`
  * Any element from any `asymp_cat` to `"V"` (if `"V"` is present)
* Place each block into `ratemat`:
  * `vax_block_dose1` -- `unvax` categories to `vaxdose_1` categories
  * `vax_block_dose2` -- `vaxprotect1` categories to `vaxdose_2` categories



# Understanding `make_vaxrate` Under Two-Dose Model 

* `vax_rate` is a vector of rates, each of which will be used to fill many rate matrix elements
* `asymp_unvax_regex` -- `"^(S|E|Ia|Ip|R)_.*unvax"`
* `asymp_vaxprotect1_regex` -- `"^(S|E|Ia|Ip|R)_.*vaxprotect1"`
  * `dose1` = `vax_prop_first_dose` * `vax_doses_per_day` / `asymp_unvax_N`
  * `dose2` = (1 - `vax_prop_second_dose`) * `vax_doses_per_day` / `asymp_unvax_N`


## (aside) Understanding `condense_state` for `return_type == "named_vector"`

Important:  This is not the same as `condense.pansim`, which is where the bottleneck _could_ be -- so the tidyverse stuff here is unlikely to be the cause of my problem.  On the other hand, this stuff could be a big contributor to the vax slowness and I bet this could be sped up simply with R refactoring (no C++ necessary)

* Pass a state vector or subset
* Seems like we are looking for states that start with these characters:
  * `states` -- `c("S", "E", "I", "H", "ICU", "R", "D")`
  * `states_regex` -- `"^S"   "^E"   "^I"   "^H"   "^ICU" "^R"   "^D"`
* Convert input to a one-row data frame (but looks like it can have more rows in some cases -- oh ... I see ... probably these cases are when we have a time-series of simulation values for each state variable) and then go tidyverse on it
  * Select using `states_regex`
  * Add an `obs_number` column -- which probably means 'time-step' in the case of condensing simulations
  * Pivot longer so that `obs_number` and `var` refer to each value
  * Split up the variable name into `state` and `subcat` by splitting on `_`
  * Make `subcat` a factor (why?)
  * Group by `obs_number` and `subcat` and summarise, which has the effect of one row per `subcat` (e.g. `vaxprotect1`), and sum up the values
  * Pivot wider so that there is one row per `obs_number` and a column for each `subcat`
  * Finally get rid of `obs_number` so that the only columns represent `subcat`
* Return the requested return type
* Answer
  * `asymp_unvax_N = sum('S_unvax', 'E_unvax', 'Ia_unvax', 'Ip_unvax', 'R_unvax')`
  * `asymp_vaxprotect1_N = sum('S_vaxprotect1', 'E_vaxprotect1', 'Ia_vaxprotect1', 'Ip_vaxprotect1', 'R_vaxprotect1')`
## 2020 April 13

* the former `aggregate` method has been split into three parts: 
    * `condense()` (collapse state variables - not called "collapse" because of namespace conflicts with `dplyr`)
	* `aggregate` (temporal aggregation)
	* `pivot` (convert from wide to long)
* calibration	
    *  moved old `calibrate` (calibration based on estimated intercept and slope from GLMM or log-linear regression) to `calibrate_slopeint`
	* new general-purpose MLE `calibrate` function (fits an arbitrary set of time series [H, ICU, D, report]; allows any parameters plus time-varying params and distributional params, link functions)
* new `forecast_ensemble()` function; uses arguments stored with `calibrate` output	
* `invlink_trans()` is now recursive
* `update.params_pansim` works for named vectors as well as individual elements


## 2020 April 15

* switched to `anytime::anydate()`: gradually moving all dates over to %Y-%m-%d or %Y-%b-%d
* EpiEstim fits
* use daily deaths rather than total (cumulative) deaths for calibration
* fit full calibration to all data, not just recent data (early case reports are needed to calibrate growth rate before the 17 March breakpoint
* alternatively try hospital-only calibration, but without the 17 March breakpoint
* major reorganization

## 2020 April 16

* continued reorganization
* process error
* improve hospital-only calibration

## 2020 April 27

* per-variable observation error
* switched to using `mle2` internally

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


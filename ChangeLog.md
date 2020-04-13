## 2020 April 13

* the former `aggregate` method has been split into three parts: 
    * `condense()` (collapse state variables - not called "collapse" because of namespace conflicts with `dplyr`)
	* `aggregate` (temporal aggregation)
	* `pivot` (convert from wide to long)
*  moved old `calibrate` (calibration based on estimated intercept and slope from GLMM or log-linear regression) to `calibrate_slopeint`; new general-purpose MLE `calibrate` method
	

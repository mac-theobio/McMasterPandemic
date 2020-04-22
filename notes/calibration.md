An attempt to get some words down to describe the current state of the calibration procedure.

## 22 April 2020

* We have abandoned the "intercept & slope" method we were using earlier (i.e. compute intercept and slope from a NB fit to exponential-phase time series, calibrate parameters to desired $r$ (= slope) and $x(0)$ by projecting backward and using the dominant eigenvector to guess relative magnitudes).  Two reasons: (1) this only works simply if we're still in exponential phase, otherwise we go back down the rabbit hole of guessing what time period to use in the statistical model; (2) there are various issues (numerical instability, rounding, etc.) in the back-projection (i.e. numerical instability can throw the initial phase off from the expected exponential growth).  We could find ways to deal with these problems, but decided to move on.
* Pure-MLE approach; calibrate mean generation interval, then pick a set of parameters (including beta0 and E0), starting time, which observables to use, etc.
    * in simplest sim examples (`run_caltest_nobreak`) no breakpoint, fit to reports only) we get bias for (at least) three reasons. 
	    * observation noise is in the wrong place in the model (should occur after convolutions, not before)
		* calibration sims start (currently) 15 days before observed time series
		* rounding
	* these *should* be the only differences between the initial sim and the calibration sims
	* if we want to fix/minimize these differences:
		* fix obs noise model
		* try fitting to hospitalization or some other response that doesn't have the convolution
		* cut the calibration data to a subset (not starting as far back)
		* use a smaller time offset (change in `make_state()` to use dominant eigenvector by default should help with this)
	* next-most-complicated sim examples (`run_caltest`)
* fits to 
    * many choices: which time series? hosp, ICU, reports, deaths?
	* which breakpoints? 
	* which auxiliary params to try to calibrate?
	* should we try to stabilize?
* known problems with calibration
    * H mismatch in model (H in model=acute care only, should include ICU as well)
	* no I -> death flow (yet)
	* nbdisp not allowed to vary across series (might need to allow Poisson?)

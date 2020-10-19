# PHAC handover

We call our model `McMasterPandemic` or "MacPan".

## Basic structure

As a quick reminder/executive summary, the basic structure is similar to many other COVID19 models

* epidemiological compartments: susceptible, exposed, presymptomatic, asymptomatic, mild and severely symptomatic
* health-care utilization: hospitalized (acute care), ICU
* tracks daily deaths, hospitalizations, ICU, case reports
* case reports are based either on a delay kernel (convolution) or a mechanistic model of testing 
* noise in observations (i.e. negative binomial variation in observation)
* all parameters can be set to vary at discrete changepoints
* additionally, Rt can be set to vary 
    * according to 'phenomenological heterogeneity' (modeling the effect of heterogeneity on population-level immunity)
	* according to a smooth time-varying curve (spline)
	* on the basis of mobility data (e.g. Google Mobility)

## Simulation

* if all parameters are specified, the model can be used for simulation/scenario projection. The 'changepoint' facility for parameters will be particularly useful; arbitrary scenarios of future changes in epidemiological parameters can easily be specified

### Poorly tested/incomplete

* mechanistic testing is fully implemented, but calibration and interpretation of 'weighting parameters' (i.e. determining relative probabilities of sampling uninfected vs presymptomatic vs asymptomatic vs symptomatic) needs work
* basic machinery for an age-structured version is present, but poorly tested
* no LTC compartment; deaths are assumed to occur only from the ICU
* no explicit vaccination model

## Calibration

By this we mean using data (on hospitalizations, deaths, case reports, and possibly testing rates or mobility data) to estimate a subset of the parameters, which can then be used for assessment (e.g. of Rt patterns) or, assuming parameters stay constant in the future, for forecasting.

This works well and has been tested allowing the baseline contact/transmission rate to vary over time according to a spline curve, using a delay model of testing or the mechanistic model of testing.

It will probably not be feasible, e.g., for age-structured models.




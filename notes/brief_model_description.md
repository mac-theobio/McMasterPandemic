
Our McMasterPandemic R package is open-source and available here:

  https://github.com/bbolker/McMasterPandemic

Our code provides tools for simulating and forecasting infectious
disease outbreaks, using compartmental epidemic models. The primary
mechanistic framework is a susceptible-exposed-infectious-removed
(SEIR) model, with additional compartments for individuals in acute
and intensive care units in hospitals. Infections can be asymptomatic,
mildly/moderately symptomatic or severely symptomatic. All symptomatic
individuals are presumed to have had a period of pre-symptomatic
infectiousness. Recovered individuals are assumed to be immune.

A fixed (fitted) proportion of deaths are assumed to occur in
hospital.

We allow for abrupt changes in transmission rate on specific dates
(associated with policy changes). The degrees of change in
transmission rate on these dates are fitted parameters.

We are able to calibrate the model to any number of observed time
series. Typically we fit simultaneously to reported cases, deaths,
hospital occupancy and ICU occupancy.

Unless otherwise stated, forecasts are based on the status quo being
maintained, i.e., we assume all parameters retain their values on the
final date of observed data.

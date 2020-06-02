## notes on model testing

### calibration

* we want to make sure that our confidence intervals have proper *coverage*, i.e. that the 95% confidence interval contains true values 95% of the time (or 90% CIs contain 90% of the values, etc.). This can only be done by (tedious) simulation, i.e. repeat the process of simulating, estimating parameters, and constructing parameter values 200x; for each forecast or parameter, count the number of times the confidence intervals contain the true value (for simplicity you can use the same true values every time, although an even more thorough analysis might

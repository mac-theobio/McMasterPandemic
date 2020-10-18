## Pros

-	sensible simulations (forecast) if we provide parameters
	- testify
	- age-structure
	- spline 

- Splines work reasonably well to real data

## Cons

- Beginning is unstable, consistently over estimating beta0
	- Start_date is very sensitive
	- testify start_date vs non_testify start_date
	- DE suggest first 10 cases or ma first 10 cases

- DE.optim Sigma vs mle2$Sigma

- What parameters can we leave out of pop_pred_samp?

- Not optimizing dispersion parameters but including them as fixed and flexibly adjusting in forecast/ensembles.





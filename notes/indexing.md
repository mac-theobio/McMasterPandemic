# Notes

* Is `pmat` a parameter that can be optimized?
* What is a parallel accumulator?

* Points to make
  * The indices do not change throughout the simulation, at least in the simple case
  * For the next step -- avoid all branches for which any of the following evaluating to `TRUE`:
    * `has_age(params)`
    * `has_testing(state = state)`
    * `has_vax(params)`
    * `has_zeta(params)` -- although I am not actually sure if this will always be off
  * If we stick to these simple cases for now, the parameters and states both fit into vectors (not matrices) -- once we figure this out, then we can consider how to unpack the matrices and pass them to C++
  * `do_step` is now entirely in C++, so I do not understand 
* Not every parameter gets optimized! How to deal with this on TMB side? -- I think this means that we have to pass in such parameters through a DATA macro -- see `opt_pars` on page 10 of the vignette
* Will `pfun` be enough for now? Check to see where/how `pfun` is used in MacPan currently -- No, it will not be enough as we also need to compute the functions of the parameters
* How to handle matrix and vector valued parameters (e.g. `pmat` for age structure?)



# Email from Mike

Hi Steve,

When I am running MacPan for PHAC, the parameters I am estimating are just beta0 and its relative values (i.e. piecewise time-varying relative scales at different time points where a policy was in place or release). Right now I am only calibrating to daily reported cases, so estimating betas are sufficient for now. 

I only have vaccination in the model and NO age structure or testing status. Irena developed the stratification, but she doesn't really run them.

I don't know how much effort it is to optimize Irena's stratification code, because it was not written in an optimized way. It was written in a convenient way for her and for her to debug. 

I would recommend sitting down and draw out the design/architect of MacPan before heading in different rabbit holes to repatch with better code. Also making this more modular would help. Let me know how I can help and always happy to sit in these meetings.

Mike 

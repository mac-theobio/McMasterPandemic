What do we need to do to reduce MacPan

- set hosp to zero (100% of cases are 'moderate')
- no asymptomatic cases
- set pre-symptomatic infectiousness equal to symptomatic infectiousness

at this point you will have a SEI^2R model ...
you can't quite make this identical to SIR because you need
the *rates* of E -> I and I1 -> I2 to become *large* which screws
things up ... it might be easier to implement an SEI^2R model in
TMB!

turn off any cleverness about using the eigenvector to establish
initial conditions, just hand it an explicit starting vector

Ali worked on this ... where is it??

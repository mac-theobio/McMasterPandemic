# Create Indices for Accessing Elements of the Rate Matrix

## Rate Matrix Models

### Warm Up Rate Matrix Model

Allow each element of the rate matrix to depend on products of
parameters (or their complements).

Let *M* = \[*M*<sub>*i*, *j*</sub>\] be the rate matrix and
*θ* = \[*θ*<sub>*k*</sub>\] be the parameter vector. The elements of *M*
are functions of the elements of *θ*. The dependence of *M* on *θ* can
only take certain forms. In particular, each element of *M* is a product
of zero or more elements of *θ* or their complements (i.e. 1 − *θ*).
Here are some examples:

-   *M*<sub>4, 5</sub> = *θ*<sub>1</sub>*θ*<sub>2</sub>
-   *M*<sub>1, 2</sub> = (1−*θ*<sub>3</sub>)*θ*<sub>4</sub>*θ*<sub>5</sub>(1−*θ*<sub>10</sub>)
-   *M*<sub>10, 3</sub> = (1−*θ*<sub>4</sub>)

In general,

*M*<sub>*i*, *j*</sub> = ∏<sub>*k*</sub>*θ*<sub>*k*</sub><sup>*x*<sub>*i**j**k*</sub></sup>(1−*θ*<sub>*k*</sub>)<sup>*y*<sub>*i**j**k*</sub></sup>

where *x*<sub>*i**j**k*</sub> (or *y*<sub>*i**j**k*</sub>) is one if
*M*<sub>*i*, *j*</sub> depends on *θ*<sub>*k*</sub> (or
1 − *θ*<sub>*k*</sub>) and zero otherwise. Note that these *x* and *y*
numbers define the dependence of *M* on *θ*, and are constants that do
not change throughout a simulation.

This model allows one to compute simple rate matrix structure where the
elements of the rate matrix are simply parameters. An example in MacPan
is the following.

    afun("Ia", "R", gamma_a)

This model also allows products of parameters and complements of
parameters, such as this example from MacPan.

    afun("Is", "ICUs", (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * gamma_s)

### More General Rate Matrix Model

#### Motivation for generalization

The above warm up model is too restrictive for some rate matrix
elements, with the most important example being the force of infection.
Here is a definition of the force of infection that is equivalent (I
think) to what MacPan uses (at least in some cases).

    afun("S", "E", sum(state[c("Ia", "Ip", "Im", "Is")] 
      * beta0 
      * c(Ca, Cp, (1 - iso_m) * Cm, (1 - iso_s) * Cs) 
      / N))

To accommodate such force of infection updates and others, we define a
more general model that allows the following additional operations.

1.  state variables can be used, in addition to parameters
2.  parameters and state variables can be inverted, in addition to
    complements
3.  products of parameters, state variables, and their inverses and
    complements can be added together

#### User Interface

Users can define the structure of a rate matrix with a list of
expressions, each determining the parameter and state dependence of a
non-zero rate matrix element. For example,

    list(
      # recovery
      list(from = "Ia", to = "R",    
           formula = ~ (gamma_a)),
      # hospitalizations
      list(from = "Is", to = "ICUs", 
           formula = ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s)),
      # force of infection
      list(from = "S",  to = "E",
           formula = ~ 
             (Ia) * (beta0) * (1/N) * (Ca) + 
             (Ip) * (beta0) * (1/N) * (Cp) + 
             (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) + 
             (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m)),
      list('etc...')
    )

The formulas allow the following operations:

-   Any element, *x*, of either the parameter or state vector can be
    placed *in parentheses* to produce a *factor* in one of the
    following three forms
    -   Identity: `(x)`
    -   Complement: `(1-x)`
    -   Inverse: `(1/x)`
-   Any number of factors can be multiplied together using `*` to
    produce a *product*
-   Any number of factors and products can be added together using `+`

#### Example

The `utilities.R` file in the folder for this document contains some
prototype functions for this approach. Using them we illustrate how rate
matrix structure can be flexibly described by the user.

    source('utilities.R')

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ratemat_struct = mk_ratemat_struct(
      rate("E", "Ia", ~ (alpha) * (sigma)),
      rate("E", "Ip", ~ (1 - alpha) * (sigma)),
      rate("Ia", "R", ~ (gamma_a)),
      rate("Ip", "Im", ~ (mu) * (gamma_p)),
      rate("Ip", "Is", ~ (1 - mu) * (gamma_p)),
      rate("Im", "R", ~ (gamma_m)),
      rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s)),
      rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s)),
      rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s)),
      rate("Is", "D", ~ (nonhosp_mort) * (gamma_s)),
      rate("ICUs", "H2", ~ (psi1)), ## ICU to post-ICU acute care
      rate("ICUd", "D", ~ (psi2)), ## ICU to death
      rate("H2", "R", ~ (psi3)), ## post-ICU to discharge
      ## H now means 'acute care' only; all H survive & are discharged
      #list(from = "H", to = "D", formula = ~ 0),
      rate("H", "R", ~ (rho)), ## all acute-care survive
      rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s)), ## assuming that hosp admissions mean *all* (acute-care + ICU)
      # force of infection
      rate("S",  "E", ~
             (Ia) * (beta0) * (1/N) * (Ca) +
             (Ip) * (beta0) * (1/N) * (Cp) +
             (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
             (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
    )

    # sanity check:
    all.equal(
      c(do_step2(state, M, params, ratemat_struct)),
      c(do_step(state, params, M, do_hazard = FALSE, do_exponential = FALSE))
    )

    ## [1] TRUE

We can parse this structure to produce information about how to update
the elements of the rate matrix. First, all of the *factors* over all
rate matrix entries in this model can be displayed.

    attributes(ratemat_struct)$all_factors

    ##             var compl invrs state_param_index factor_index
    ## 1         alpha FALSE FALSE                20            1
    ## 2         sigma FALSE FALSE                21            2
    ## 3         alpha  TRUE FALSE                20            3
    ## 4       gamma_a FALSE FALSE                22            4
    ## 5            mu FALSE FALSE                28            5
    ## 6       gamma_p FALSE FALSE                25            6
    ## 7            mu  TRUE FALSE                28            7
    ## 8       gamma_m FALSE FALSE                23            8
    ## 9  nonhosp_mort  TRUE FALSE                31            9
    ## 10         phi1 FALSE FALSE                34           10
    ## 11      gamma_s FALSE FALSE                24           11
    ## 12         phi1  TRUE FALSE                34           12
    ## 13         phi2  TRUE FALSE                35           13
    ## 14         phi2 FALSE FALSE                35           14
    ## 15 nonhosp_mort FALSE FALSE                31           15
    ## 16         psi1 FALSE FALSE                36           16
    ## 17         psi2 FALSE FALSE                37           17
    ## 18         psi3 FALSE FALSE                38           18
    ## 19          rho FALSE FALSE                26           19
    ## 20           Ia FALSE FALSE                 3           20
    ## 21        beta0 FALSE FALSE                15           21
    ## 22            N FALSE  TRUE                29           22
    ## 23           Ca FALSE FALSE                16           23
    ## 24           Ip FALSE FALSE                 4           24
    ## 25           Cp FALSE FALSE                17           25
    ## 26           Im FALSE FALSE                 5           26
    ## 27           Cm FALSE FALSE                18           27
    ## 28        iso_m  TRUE FALSE                32           28
    ## 29           Is FALSE FALSE                 6           29
    ## 30           Cs FALSE FALSE                19           30

Here is an example where this input,

    ratemat_struct[[5]][c('from', 'to', 'formula')]

    ## $from
    ## [1] "Ip"
    ## 
    ## $to
    ## [1] "Is"
    ## 
    ## $formula
    ## ~(1 - mu) * (gamma_p)

parses as this,

    ratemat_struct[[5]][c('factors')]

    ## $factors
    ##   product_index     var compl invrs state_param_index factor_index
    ## 1             1      mu  TRUE FALSE                28            7
    ## 2             1 gamma_p FALSE FALSE                25            6

#### Mathematical Model Description

We make a series of transformations from the state vector, *s*, and
parameter vector to the rate matrix, *M*. Given how TMB works, we need
to distinguish parameters that will be optimized, *θ*, from those that
will remain constant, *η* (or those that are derived from parameters
being optimized, e.g. time-varying??, still unsure how time varying
machinery is working with run\_sim\_break or will work with
run\_sim\_loglin). We define two intermediate vectors as well: the
factor vector, *v*, and the product vector, *u*. At a high level we take
the elements of the state vector and two parameter vectors into the
non-zero elements of the rate matrix as follows.
{*s*, *η*, *θ*} → *v* → *u* → *M*

The factor vector, *v*, is a function of *s* and *θ*. The dependence is
simple in that each scalar element of the factor vector can be one of
the following.

-   \[identity\] an element, *x*, of either *s* or *θ*
-   \[complement\] the complement of this element, 1 − *x*
-   \[inverse\] the inverse of this element, 1/*x*

In the future we can generalize this by adding more operations, but for
now identity, complement, and inverse should get us pretty far –
although it will make age structure awkward given that such problems are
more naturally handled with matrix operations like Kronecker products
and sweeps.

The product vector, *u*, is a function of *v*. The dependence is simple
in that each element of the product vector is the product of one or more
elements of the factor vector, *v*.

The non-zero elements of the rate matrix, *M*, can be any one of the
following.

-   An element of the factor vector, *v*
-   An element of the product vector, *u*
-   The sum of one or more elements in the product vector, *u*

Explicitly expressing the force of infection in terms of this model we
have.

*s* = \[*I*<sub>*a*</sub>,*I*<sub>*p*</sub>,*I*<sub>*m*</sub>,*I*<sub>*s*</sub>\]

*θ* = \[*β*<sub>0</sub>,*C*<sub>*a*</sub>,*C*<sub>*p*</sub>,*C*<sub>*m*</sub>,*C*<sub>*s*</sub>,*i**s**o*<sub>*m*</sub>,*i**s**o*<sub>*s*</sub>,*N*\]

*u* = \[*s*<sub>1</sub>,*s*<sub>2</sub>,*s*<sub>3</sub>,*s*<sub>4</sub>,*θ*<sub>1</sub>,*θ*<sub>2</sub>,*θ*<sub>3</sub>,*θ*<sub>4</sub>,*θ*<sub>5</sub>,1−*θ*<sub>6</sub>,1−*θ*<sub>7</sub>,1/*θ*<sub>8</sub>\]

*v* = \[*u*<sub>1</sub>*u*<sub>5</sub>*u*<sub>6</sub>*u*<sub>12</sub>,*u*<sub>2</sub>*u*<sub>5</sub>*u*<sub>7</sub>*u*<sub>12</sub>,*u*<sub>3</sub>*u*<sub>5</sub>*u*<sub>10</sub>*u*<sub>8</sub>*u*<sub>12</sub>,*u*<sub>4</sub>*u*<sub>5</sub>*u*<sub>11</sub>*u*<sub>9</sub>*u*<sub>12</sub>\]

And the non-zero element of *M* determining the rate of flow from S to E
is the sum of the elements in the products vector, *v*. Typically this
summation will be taken over a subset of the elements of the products
vector, but we simplified this example to include on elements that are
involved in the force of infection computation.

## Roadmap

### Avoid Copying on the R-Side

Although the overall motivation for this work is to facilitate C++/TMB
refactoring, it would be good to have the R-side version as efficient as
possible.

### Constants

Could be convenient to just allow users to hard code constant rate
matrix entries in the rate matrix structure. This is an example from
MacPan.

    afun("H", "D", 0)

I’m surprised this is necessary though because I have been assuming that
unspecified rate matrix elements imply zero.

### Regex Matching

Not thinking about this at all, but probably should. For example, could
`pfun` return multiple rate matrix positions? If so that would screw up
the current index machinery.

### Index Permutations

Permutation of the parameter, state, and factor vectors for
computational efficiency. We don’t want to have to sum together products
that depend on elements of the state vector that are far apart from each
other in memory. This is similar to permutations of sparse matrices in
the Matrix package and Eigen.

### Model Specification in Terms of Matrices and Vectors

The model above treats all parameters and state variables as scalars.
But in many cases it is more natural to consider vector- and
matrix-valued parameters and states. For example, the contact matrix
(`pmat`) that is used by MacPan in models of age-structure.

As far as I can tell models such as age structure *could* be defined in
terms of the model above (because matrix operations are composed of
products and sums of products), but it would just require large numbers
of tedious entries. There are two ways around this tedium. First, we
could add matrix operations to the list of allowable operations in the
rate matrix update model. Second, we could create convenience
model-specification utilities that would automatically expand to the
notation used in the above model. I’m not sure which one I like better,
but the first one could take advantage of sparse matrix optimizations
from Eigen etc. On the other hand, maybe we can get similar performance
with smart index permutations on the R-side?

## Unstructured Notes

When coding this up, we do not need to explicitly store all of the zeros
in *x* and *y*. Instead we just track what elements of *M* depend on
what elements of *θ* and whether/how the dependence involves complements
of the values of *θ*.

Before defining the interface, assume that we have a rate matrix on the
R-side with named rows and columns, and a named parameter vector also on
the R-side with named elements.

    theta <- McMasterPandemic:::read_params("ICU1.csv")
    state = McMasterPandemic:::make_state(params=theta)
    M <- McMasterPandemic:::make_ratemat(params=theta, state=state, sparse=TRUE)
    print(theta)

    ##        beta0           Ca           Cp           Cm           Cs        alpha 
    ## 1.000000e+00 6.666667e-01 1.000000e+00 1.000000e+00 1.000000e+00 3.333333e-01 
    ##        sigma      gamma_a      gamma_m      gamma_s      gamma_p          rho 
    ## 1.923077e-01 1.428571e-01 1.428571e-01 1.748252e-01 2.000000e+00 1.000000e-01 
    ##        delta           mu            N           E0 nonhosp_mort        iso_m 
    ## 0.000000e+00 9.560000e-01 1.000000e+06 5.000000e+00 0.000000e+00 0.000000e+00 
    ##        iso_s         phi1         phi2         psi1         psi2         psi3 
    ## 0.000000e+00 7.600000e-01 5.000000e-01 5.000000e-02 1.250000e-01 2.000000e-01 
    ##       c_prop c_delay_mean   c_delay_cv    proc_disp         zeta 
    ## 1.000000e-01 1.100000e+01 2.500000e-01 0.000000e+00 0.000000e+00

    print(M)

    ## 14 x 14 sparse Matrix of class "dgTMatrix"

    ##    [[ suppressing 14 column names 'S', 'E', 'Ia' ... ]]

    ##                                                                               
    ## S    . 1.666667e-06 .          .         .     .     .         .    .         
    ## E    . .            0.06410256 0.1282051 .     .     .         .    .         
    ## Ia   . .            .          .         .     .     .         .    .         
    ## Ip   . .            .          .         1.912 0.088 .         .    .         
    ## Im   . .            .          .         .     .     .         .    .         
    ## Is   . .            .          .         .     .     0.1328671 .    0.02097902
    ## H    . .            .          .         .     .     .         .    .         
    ## H2   . .            .          .         .     .     .         .    .         
    ## ICUs . .            .          .         .     .     .         0.05 .         
    ## ICUd . .            .          .         .     .     .         .    .         
    ## D    . .            .          .         .     .     .         .    .         
    ## R    . .            .          .         .     .     .         .    .         
    ## X    . .            .          .         .     .     .         .    .         
    ## V    . .            .          .         .     .     .         .    .         
    ##                                            
    ## S    .          .     .         .         .
    ## E    .          .     .         .         .
    ## Ia   .          .     0.1428571 .         .
    ## Ip   .          .     .         .         .
    ## Im   .          .     0.1428571 .         .
    ## Is   0.02097902 .     .         0.1328671 .
    ## H    .          .     0.1000000 .         .
    ## H2   .          .     0.2000000 .         .
    ## ICUs .          .     .         .         .
    ## ICUd .          0.125 .         .         .
    ## D    .          .     .         .         .
    ## R    .          .     .         .         .
    ## X    .          .     .         .         .
    ## V    .          .     .         .         .

    state_names = row.names(M)
    par_names = names(theta)
    print(state_names)

    ##  [1] "S"    "E"    "Ia"   "Ip"   "Im"   "Is"   "H"    "H2"   "ICUs" "ICUd"
    ## [11] "D"    "R"    "X"    "V"

    print(par_names)

    ##  [1] "beta0"        "Ca"           "Cp"           "Cm"           "Cs"          
    ##  [6] "alpha"        "sigma"        "gamma_a"      "gamma_m"      "gamma_s"     
    ## [11] "gamma_p"      "rho"          "delta"        "mu"           "N"           
    ## [16] "E0"           "nonhosp_mort" "iso_m"        "iso_s"        "phi1"        
    ## [21] "phi2"         "psi1"         "psi2"         "psi3"         "c_prop"      
    ## [26] "c_delay_mean" "c_delay_cv"   "proc_disp"    "zeta"

    pfun <- McMasterPandemic:::pfun

    make_beta_light <- function(state, params, full = TRUE) {
        # does not use 'full' (argument)
        # does not use 'testcats' (below)
        Icats <- c("Ia", "Ip", "Im", "Is")
        testcats <- c("_u", "_p", "_n", "_t")
        Icat_prop_vec <- with(
            as.list(params),
            c(Ca, Cp, (1 - iso_m) * Cm, (1 - iso_s) * Cs)
        )
        names(Icat_prop_vec) <- Icats
        beta_0 <- with(as.list(params), beta0 * Icat_prop_vec / N)
        beta <- setNames(numeric(length(state)), names(state))
        beta[names(beta_0)] <- beta_0
        return(beta)
    }

    update_ratemat_light <- function(ratemat, state, params, testwt_scale = "N") {
      # testwt_scale not needed
      ratemat[pfun("S", "E", ratemat)] <- update_foi(state, params, make_beta(state, params))
    }

    update_foi_light <- function(state, params, beta) {
      # params not needed
      foi <- sum(state * beta)
    }

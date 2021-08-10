# Dependence of State Transition Rates on Parameters and State Variables

Currently MacPan developers specify new models by writing R code that
transforms parameters and state variables into state transition rates
(`make_ratemat`, `update_ratemat`). As we build more of MacPan in `C++`
we would like to maintain the ability of epidemiologically-focused
developers to define new models in R as much as possible. At the same
time we would like to have `C++` developers focused on computational
efficiency rather than details of epidemiological models. In order to
achieve this separation of concerns, we will need to formalize a spec
for a general model of how state transition rates depend on parameters
and state variables. This spec will make it easier for the two types of
developers to work together. We will version this spec so that we can
build it in steps, making sure that we don’t over-engineer it.

This is a working document that provides background on current thinking,
spec definitions (currently only drafts), and a roadmap for planning and
remembering details that are not yet sorted out. Over time items from
the roadmap will move into defined spec versions. To update it use
`make` from this folder.

## Versioned Rate Matrix Models

There are currently no versioned rate matrix model specs. See below for
experimental models.

## Experimental Rate Matrix Models

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

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    rs = mk_ratemat_struct(
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
      c(do_step2(state, M, params, rs)),
      c(do_step(state, params, M, do_hazard = FALSE, do_exponential = FALSE))
    )

    ## [1] TRUE

We can parse this structure to produce information about how to update
the elements of the rate matrix. First, all of the *factors* over all
rate matrix entries in this model can be displayed.

    attributes(rs)$all_factors

    ##             var compl invrs sim_updates opt_updates state_param_index
    ## 1         alpha FALSE FALSE       FALSE       FALSE                20
    ## 2         sigma FALSE FALSE       FALSE       FALSE                21
    ## 3         alpha  TRUE FALSE       FALSE       FALSE                20
    ## 4       gamma_a FALSE FALSE       FALSE       FALSE                22
    ## 5            mu FALSE FALSE       FALSE       FALSE                28
    ## 6       gamma_p FALSE FALSE       FALSE       FALSE                25
    ## 7            mu  TRUE FALSE       FALSE       FALSE                28
    ## 8       gamma_m FALSE FALSE       FALSE       FALSE                23
    ## 9  nonhosp_mort  TRUE FALSE       FALSE       FALSE                31
    ## 10         phi1 FALSE FALSE       FALSE       FALSE                34
    ## 11      gamma_s FALSE FALSE       FALSE       FALSE                24
    ## 12         phi1  TRUE FALSE       FALSE       FALSE                34
    ## 13         phi2  TRUE FALSE       FALSE       FALSE                35
    ## 14         phi2 FALSE FALSE       FALSE       FALSE                35
    ## 15 nonhosp_mort FALSE FALSE       FALSE       FALSE                31
    ## 16         psi1 FALSE FALSE       FALSE       FALSE                36
    ## 17         psi2 FALSE FALSE       FALSE       FALSE                37
    ## 18         psi3 FALSE FALSE       FALSE       FALSE                38
    ## 19          rho FALSE FALSE       FALSE       FALSE                26
    ## 20           Ia FALSE FALSE        TRUE        TRUE                 3
    ## 21        beta0 FALSE FALSE       FALSE        TRUE                15
    ## 22            N FALSE  TRUE       FALSE       FALSE                29
    ## 23           Ca FALSE FALSE       FALSE       FALSE                16
    ## 24           Ip FALSE FALSE        TRUE        TRUE                 4
    ## 25           Cp FALSE FALSE       FALSE       FALSE                17
    ## 26           Im FALSE FALSE        TRUE        TRUE                 5
    ## 27           Cm FALSE FALSE       FALSE       FALSE                18
    ## 28        iso_m  TRUE FALSE       FALSE       FALSE                32
    ## 29           Is FALSE FALSE        TRUE        TRUE                 6
    ## 30           Cs FALSE FALSE       FALSE       FALSE                19
    ##    factor_index
    ## 1             1
    ## 2             2
    ## 3             3
    ## 4             4
    ## 5             5
    ## 6             6
    ## 7             7
    ## 8             8
    ## 9             9
    ## 10           10
    ## 11           11
    ## 12           12
    ## 13           13
    ## 14           14
    ## 15           15
    ## 16           16
    ## 17           17
    ## 18           18
    ## 19           19
    ## 20           20
    ## 21           21
    ## 22           22
    ## 23           23
    ## 24           24
    ## 25           25
    ## 26           26
    ## 27           27
    ## 28           28
    ## 29           29
    ## 30           30

Here is an example where this input,

    rs$Ip_to_Is[c('from', 'to', 'formula')]

    ## $from
    ## [1] "Ip"
    ## 
    ## $to
    ## [1] "Is"
    ## 
    ## $formula
    ## ~(1 - mu) * (gamma_p)

parses as this,

    rs$Ip_to_Is$factors

    ##   product_index     var compl invrs sim_updates.x opt_updates.x sim_updates.y
    ## 1             1      mu  TRUE FALSE         FALSE         FALSE         FALSE
    ## 2             1 gamma_p FALSE FALSE         FALSE         FALSE         FALSE
    ##   opt_updates.y state_param_index factor_index
    ## 1         FALSE                28            7
    ## 2         FALSE                25            6

More generally, all elements of a `ratemat-struct` object are
`rate-struct` objects and contain the following elements:

    names(rs$Is_to_H)

    ## [1] "from"            "to"              "formula"         "factors"        
    ## [5] "ratemat_indices"

The `from` and `to` elements give the name of the associated states,
`formula` gives the expression determining the rate, `factors` gives a
data frame where each row is a *factor* (as defined above) required to
compute the rate, and `ratemat_indices` gives the row and column index
of the rate in the rate matrix.

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

The factor vector, *v*, is a function of *s*, *η*, and *θ*. The
dependence is simple in that each scalar element of the factor vector
can be one of the following.

-   \[identity\] an element, *x*, of *s*, *η*, or *θ*
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

### Parameters to optimize vs other parameters

Because of how TMB is used to produce objective functions for
optimization, the parameters to optimize/calibrate must be separated
from the other parameters. Currently I have ignored this complexity, but
I think we don’t need to think about it until we actually implement this
on the `C++` side. This separation of parameters will still need to be
done thoughtfully though, because it might have implications for
accessing memory efficiently during simulation steps.

### Index Permutations

Permutation of the parameter, state, and factor vectors for
computational efficiency. We don’t want to have to sum together products
that depend on elements of the state vector that are far apart from each
other in memory. This is similar to permutations of sparse matrices in
the Matrix package and Eigen.

Made some progress here – we can visualize the dependence of state
transition rates on variables (either parameters or state variables) by
plotting a heat map with transitions on one axis and variables on the
other. The colour of the heatmap depends on the number of times a
variable appears in the expression for the transition rate.

    m = as.matrix(rs)

    ## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

    ## `summarise()` has grouped output by 'var'. You can override using the `.groups` argument.

    ca = corresp(m)
    heatmap(m[order(ca$rscore), order(ca$cscore)], 
            Colv = NA, Rowv = NA, scale = "none")

![](indexing_files/figure-markdown_strict/unnamed-chunk-11-1.png)

I’ve ordered the axes using correspondence analysis, which is a
statistical technique for finding block structure in matrices. Block
structure might be relevant because such structure means that we can
access parameters in chunks, which could be useful for optimizing on the
`C++` side.

### Time Varying Parameters

The model above does not consider time varying parameters, but this will
need to be done obviously. My current sense is that for
break-point-style parameter changes should be handled when we separate
parameters to be optimized from the other parameters (see above).
Log-linear-style parameter changes, should probably be handled when we
introduce model specification in terms of matrices and vectors (see
below).

Now that I’ve read the MacPan manuscript I feel better about the
log-linear model of the transmission rate, *β*. A generalization to
something beyond the transmission rate is the following.

log (*θ*<sub>*i*</sub>(*t*)) = log (*θ*<sub>*i*</sub><sup>(0)</sup>) + **X**

Each row of **X** corresponds to

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

One issue with the first option is that there is only one kind of scalar
multiplication, but several kinds of matrix multiplication – so we can’t
just overload the `*` operator (although we probably should as is done
in `col_multiply`) and will need to introduce explicit matrix operators
(e.g. `%*%`, `%x%`). An issue with Kronecker-type operations (as are
used in age-structured models) is that they often involve matrix
transposing.

There are some relevant notes on this topic here:
<https://github.com/mac-theobio/McMasterPandemic/blob/master/notes/structure.md#general-concerns>

## Steve’s Notes on Meeting with Ben (2021-08-09)

-   very important to be able to specify multiple rates in the model
    specification machinery – this is what the regex stuff is doing now
-   Additional motivation: figuring out indexing ahead of simulation
    time is also likely to be an efficiency win.
-   3 components of the model that involve updates of state variables
    (TODO: verify that the model works for these cases at a minimum):
    -   FOI
    -   testing – vec of weights that determines prob of being tested
        depending on state
    -   vaccination – we know realized vaccination rates, but need to
        scale to per-capita
-   Kronecker products may only be convenient for producing a matrix,
    rather than efficient at updating it
-   most elements of the rate matrix only need to be touched once
    (update/make rate matrix)
    -   this needs to go into the indexing framework
    -   by defining in terms of per-capita flows means that most
        elements can stay fixed
        1.  once at beginning of opt run
        2.  every time a new set of params are updated
        3.  done within the simulation (e.g. state-dependence and time
            varying parameters
-   logs of parameters – think more about them
-   block structure of rate dependence could have a connection with
    identifiability – there exists a literature on identifiability of
    dynamical systems models that boarders on practical
-   `TMB::MakeADFun` has a map argument to specify that parameters can
    be fixed at their starting values or set groups of parameters to be
    equal – should use this instead of switching back and forth between
    the `PARAMETER` and `DATA` macros

[![Superseded](https://img.shields.io/static/v1.svg?label=Lifecycle&message=Superseded&color=red)](https://canmod.net/misc/flex_specs)

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
spec definitions, and a roadmap for planning and remembering details
that are not yet sorted out. Over time items from the roadmap will move
into defined spec versions. To update it use `make` from this folder.

## Versioned Model Specs

Using semantic versioning <https://semver.org/>. All versions of the
spec that have not been implemented in a released McMasterPandemic R
package have pre-release `-alpha` identifiers
(<https://semver.org/#spec-item-9>).

### 0.0.1-alpha

#### Summary of Capabilities (0.0.1-alpha)

Update the elements of a rate matrix on the `TMB`/`C++` side using a
restricted set of operations. The update formula must be specified
separately for each element.

#### Assumptions (0.0.1-alpha)

-   `vaxify`, `ageify`, and `testify` are all `FALSE` for all
    `McMasterPandemic` function calls
-   No parameter vary throughout a simulation
-   A `param_pansim` and a `state_pansim` object exists or can be
    constructed
-   The results of running the parameter and state vectors against
    `make_ratemat` are stored

#### Interface (0.0.1-alpha)

Users can define the structure of a rate matrix with a list of
expressions, each determining the parameter and state dependence of a
single non-zero rate matrix element. The formulas allow the following
operations:

-   Any element, *x*, of either the parameter or state vector can be
    placed *in parentheses* to produce a *factor* in one of the
    following three forms
    -   Identity: `(x)`
    -   Complement: `(1-x)`
    -   Inverse: `(1/x)`
-   Any number of factors can be multiplied together using `*` to
    produce a *product*
-   Any number of factors and products can be added together using `+`

#### Data Structure (0.0.1-alpha)

Five integer vectors must be passed to TMB to define a rate matrix model

Each of the elements of three of these vectors correspond to each of the
non-zero elements of the rate matrix. These three vectors are defined as
follows.

-   `from`: Vector of indices pointing to the rows of the rate matrix to
    be updated

##### `from`

-   *Description*: Vector over the of 1-based indices pointing to the
    rows corresponding to the of the rate matrix to be updated. This
    index identifies a compartment in the model *from* which individuals
    flow. The `from` vector and `to` vector ‘line-up’, in the sense that
    the *i*th position in `from` indexes the state that is flowing to
    the state indexed by the *i*th position in `to`
-   *Length*: Number of scalar-valued states. If the model contains
    vector or matrix valued states, then each element of these vectors
    and matrices is counted separately
-   *Example*: 2, 2, 3, 4, 4, 5, 6, 6, 6, 6, 9, 10, 8, 7, 6, 1

##### `to`

##### `count`

##### `spi`

##### `modifier`

### In-Development Model Specs

## Background

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

    mk_ratemat_struct(
      # recovery
      rate(from = "Ia", to = "R",    
           formula = ~ (gamma_a)),
      # hospitalizations
      rate(from = "Is", to = "ICUs", 
           formula = ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s)),
      # force of infection
      rate(from = "S",  to = "E",
           formula = ~ 
             (Ia) * (beta0) * (1/N) * (Ca) + 
             (Ip) * (beta0) * (1/N) * (Cp) + 
             (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) + 
             (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m)),
      rate('etc...')
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

    ##             var compl invrs sim_updt opt_updt var_indx factor_indx
    ## 1         alpha FALSE FALSE    FALSE    FALSE       20           1
    ## 2         sigma FALSE FALSE    FALSE    FALSE       21           2
    ## 3         alpha  TRUE FALSE    FALSE    FALSE       20           3
    ## 4       gamma_a FALSE FALSE    FALSE    FALSE       22           4
    ## 5            mu FALSE FALSE    FALSE    FALSE       28           5
    ## 6       gamma_p FALSE FALSE    FALSE    FALSE       25           6
    ## 7            mu  TRUE FALSE    FALSE    FALSE       28           7
    ## 8       gamma_m FALSE FALSE    FALSE    FALSE       23           8
    ## 9  nonhosp_mort  TRUE FALSE    FALSE    FALSE       31           9
    ## 10         phi1 FALSE FALSE    FALSE    FALSE       34          10
    ## 11      gamma_s FALSE FALSE    FALSE    FALSE       24          11
    ## 12         phi1  TRUE FALSE    FALSE    FALSE       34          12
    ## 13         phi2  TRUE FALSE    FALSE    FALSE       35          13
    ## 14         phi2 FALSE FALSE    FALSE    FALSE       35          14
    ## 15 nonhosp_mort FALSE FALSE    FALSE    FALSE       31          15
    ## 16         psi1 FALSE FALSE    FALSE    FALSE       36          16
    ## 17         psi2 FALSE FALSE    FALSE    FALSE       37          17
    ## 18         psi3 FALSE FALSE    FALSE    FALSE       38          18
    ## 19          rho FALSE FALSE    FALSE    FALSE       26          19
    ## 20           Ia FALSE FALSE     TRUE     TRUE        3          20
    ## 21        beta0 FALSE FALSE    FALSE     TRUE       15          21
    ## 22            N FALSE  TRUE    FALSE    FALSE       29          22
    ## 23           Ca FALSE FALSE    FALSE    FALSE       16          23
    ## 24           Ip FALSE FALSE     TRUE     TRUE        4          24
    ## 25           Cp FALSE FALSE    FALSE    FALSE       17          25
    ## 26           Im FALSE FALSE     TRUE     TRUE        5          26
    ## 27           Cm FALSE FALSE    FALSE    FALSE       18          27
    ## 28        iso_m  TRUE FALSE    FALSE    FALSE       32          28
    ## 29           Is FALSE FALSE     TRUE     TRUE        6          29
    ## 30           Cs FALSE FALSE    FALSE    FALSE       19          30

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

    ##       var compl invrs prod_indx sim_updt opt_updt var_indx factor_indx
    ## 1      mu  TRUE FALSE         1    FALSE    FALSE       28           7
    ## 2 gamma_p FALSE FALSE         1    FALSE    FALSE       25           6

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

We make a series of transformations from the state vector and parameter
vector to the rate matrix, *M*. We concatenate the state vector and
parameter vector into a single *variable* vector called *v*. We define
two intermediate vectors as well: the factor vector, *u*, and the
product vector, *w*. At a high level we take the elements of the
variable vector into the non-zero elements of the rate matrix as
follows.
*v* → *u* → *w* → *M*

The factor vector, *u*, is a function of the variable vector, *v*. The
dependence is simple in that each scalar element of the factor vector
can be one of the following.

-   \[identity\] an element, *x*, of *v*
-   \[complement\] the complement of this element, 1 − *x*
-   \[inverse\] the inverse of this element, 1/*x*

In the future we can generalize this by adding more operations, but for
now identity, complement, and inverse should get us pretty far –
although it will make age structure awkward given that such problems are
more naturally handled with matrix operations like Kronecker products
and sweeps.

The product vector, *w*, is a function of *u*. The dependence is simple
in that each element of the product vector is the product of one or more
elements of the factor vector, *u*.

The non-zero elements of the rate matrix, *M*, can be any one of the
following.

-   An element of the factor vector, *u*
-   An element of the product vector, *w*
-   The sum of one or more elements in the product vector, *w*

Explicitly expressing the force of infection in terms of this model we
have.

*v* = \[*I*<sub>*a*</sub>,*I*<sub>*p*</sub>,*I*<sub>*m*</sub>,*I*<sub>*s*</sub>,*β*<sub>0</sub>,*C*<sub>*a*</sub>,*C*<sub>*p*</sub>,*C*<sub>*m*</sub>,*C*<sub>*s*</sub>,*i**s**o*<sub>*m*</sub>,*i**s**o*<sub>*s*</sub>,*N*\]

*u* = \[*v*<sub>1</sub>,*v*<sub>2</sub>,*v*<sub>3</sub>,*v*<sub>4</sub>,*v*<sub>5</sub>,*v*<sub>6</sub>,*v*<sub>7</sub>,*v*<sub>8</sub>,*v*<sub>9</sub>,1−*v*<sub>10</sub>,1−*v*<sub>11</sub>,1/*v*<sub>12</sub>\]

*w* = \[*u*<sub>1</sub>*u*<sub>5</sub>*u*<sub>6</sub>*u*<sub>12</sub>,*u*<sub>2</sub>*u*<sub>5</sub>*u*<sub>7</sub>*u*<sub>12</sub>,*u*<sub>3</sub>*u*<sub>5</sub>*u*<sub>10</sub>*u*<sub>8</sub>*u*<sub>12</sub>,*u*<sub>4</sub>*u*<sub>5</sub>*u*<sub>11</sub>*u*<sub>9</sub>*u*<sub>12</sub>\]

And the non-zero element of *M* determining the rate of flow from S to E
is the sum of the elements in the product vector, *w*. Typically this
summation will be taken over a subset of the elements of the products
vector, but we simplified this example to include only elements that are
involved in the force of infection computation.

## Roadmap

### Common Factors

It should be possible to express the force of infection in this way:

    (beta0) * (1/N) * (
      (Ia) * (Ca) +
      (Ip) * (Cp) +
      (Im) * (Cm) * (1-iso_m) +
      (Is) * (Cs) * (1-iso_m))

Rather than this way:

    (Ia) * (beta0) * (1/N) * (Ca) +
    (Ip) * (beta0) * (1/N) * (Cp) +
    (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
    (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m)

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

![](/Users/stevenwalker/Development/McMasterPandemic/notes/indexing/doc-output/indexing_files/figure-markdown_strict/unnamed-chunk-12-1.png)

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

log (*θ*<sub>*i*</sub>(*t*)) = log (*θ*<sub>*i*</sub><sup>(0)</sup>) + **X****c**

Each row of **X** corresponds to … TODO

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

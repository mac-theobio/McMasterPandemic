# Create Indices for Accessing Elements of the Rate Matrix

## Dependence of the Rate Matrix on the Parameter Vector

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

## User Interface

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

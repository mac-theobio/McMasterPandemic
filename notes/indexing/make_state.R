#### Make State Vector (observed)

1.  Define an initial state vector $s$, with all zeros. In the simple
vax example this is a length-70 vector of zeros.

2.  Define an initializing state vector, $\tilde{s}$, with all zeros. In
the simple vax example this is a length-60 vector of zeros. The
difference between 70 and 60 is that $\tilde{s}$ does not include
the ten X and V epi-states.

3.  Create a copy, $\tilde{\theta}$, of the parameter vector, $\theta$,
and replace some of the values with constant user-supplied values.
In the simple vax example these replacements are.

-   `N = 1`
-   `vax_doses_per_day = 0`
-   `vax_response_rate = 0`
-   `vax_response_rate_R = 0`
-   `E_0 = 1e-5`

4.  Replace some of the state vector elements with elements of the
parameter vector. In the simple vax example these replacements are.

-   `S_unvax = 1 - E_0`
-   `E_unvax = E_0` -- but this should just be in the parameter vector
instead of requiring the complement (Notes from meeting with Ben:
                                       leave S_unvax at 1 = N, and E_unvax = E_0 where E_0 is small)

5.  Make the rate matrix from $\tilde{s}$ and $\tilde{\theta}$, using
`from`, `to`, `count`, `spi`, and `modifier` (with these index and
                                              count vectors possibly adjusted to account for the lack of X and V
                                              states). (Notes: probably don't need to drop these before
    calculating the Jacobian, but definitely before the eigenvector)

6.  Run a simulation for (by default) 100 iterations with
    $\tilde{\theta}$ and starting from $\tilde{s}$. The formulas for
    calculating the outflows must be generalized relative to the
    existing version.

7.  Let $f_i^{(out)}$ be the $i$th element of the outflow vector, which
    is computed as follows. $$
    f_i^{(out)} = \sum_{j \in \Omega_i}F_{ij}
    $$

8.  Where $\Omega_i$ is a subset of the indices into the state vector.
    Typically there are only a small number of unique index sets, with
    one $\Omega_i$ for many values of $i$ -- for this reason the data
    structure for these indices only need to track each unique index
    set.

9.  The resulting final state vector, $\tilde{s}$, is saved and
    normalized such that $\sum_i\tilde{s}_i = 1$.

```{=html}
<!-- -->
```
    ret -- eigenvector -- selected states of the eigenvector and normalized over those selected states
               E           Ia           Ip           Im           Is            H           H2
    6.632989e-01 1.066418e-01 3.769675e-02 1.807735e-01 7.702539e-03 2.875948e-03 5.794981e-05
            ICUs         ICUd
    5.283317e-04 4.242892e-04

    state -- disease-free-state
          S       E      Ia      Ip      Im      Is       H      H2    ICUs    ICUd       D       R
    0.99999 0.00001 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000

    drop_vars <- c("S", "date", "t", "D", "foi", "R", "X", "N", "P", "V")

    E_values = round(evec * E0) -- a selection
       E   Ia   Ip   Im   Is    H   H2 ICUs ICUd
       3    1    0    1    0    0    0    0    0

    S_values = N - sum(E_values)
    [1] 999995

    state
    S      E     Ia     Ip     Im     Is      H     H2   ICUs   ICUd      D      R      X      V
    999995 3      1      0      1      0      0      0      0      0      0      0      0      0

Once you get an eigenvector: - Uninfected states stay at disease-free
equilibrium - Infected states gets total number of initial infected
times the scaled eigenvector

-   `disease_free_state`

    -   vector for jac calc
    -   require that `sum(disease_free_state) == 1`

-   `is_infected_state` (should be true when `disease_free_state != 0`)

-   `initial_state[!is_infected_state] *= N - E_0`

-   `initial_state[is_infected_state] = eigenvector[is_infected_state] * E_0`

-   `initial_state = smart_round(initial_state)`

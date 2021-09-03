# Designing a data structure to be passed under TMB to calculate rate matrix in C++
 In this draft, we propose a data structure that meets the requirements of [specs from 
 0.0.1 to 0.0.3](https://canmod.net/misc/flex_specs) (perhaps even including 0.0.4) 
 written by Steven. The proposed data structure described below is both extension and 
 alteration of previous proposals. This new proposal addresses the time-varying parameters 
 in the calculation of rate matrix.
 
 ## Data Structure
The parts of the data structure that remain unchanged are:
 - *sp*: state vector 
 concatenated with parameter vector (the same as in specs of 0.0.1 and 0.0.2). A 
 time-varying parameter still holds one place in the vector rather than expanding *sp* to 
 hold multiple values of time-varying factors.
 - *from*,*to*, *count*, *spi* and 
 *modifier* are defined in the same as in previous proposals.
The parts of the data structure that are added, removed or alterred are:
 - *upateidx* (**added**): a vector of indices into vectors *from*, *to*, and *count* of those elements 
 in the rate matrix that need to be updated. It includes the indices of the elements that 
 depend on either the state vector (0.0.2), or time-varying parameters (0.0.3 or 0.0.4), 
 or both. At each simulation step, we will update only those elements specified by 
 updateidx. The drawback of this design is that elements that only depends on piece-wise 
 parameters (0.0.3) are updated at each step although they only need to be updated at 
 certain breaks.
 - *update_indices* (**removed**): a second set of *from* *to*, etc for 
 need-to-update elements proposed in spec 0.0.2 is removed in this proposal because the 
 newly introduced *upateidx* provides the information needed.
 - *breaks* (**added**): 
 vector of all the breaks.
 - *count_of_tv_at_breaks (**added**): vector of number of 
 time-varying parameters that change at each break.
 - *tv_spi* (**added**): vector of 
 indices into vector *sp* of those time-varying parameters in the order of breaks as major 
 and parameters as minor.
 - *tv_val* (**added**): vector of new values corresponding to 
 *tv_spi*.

To help you understand the proposed data structure, lets take an example. Suppose the 
time-varying parameters are *beta0* and *mu*. Their indices in *sp* are 15 and 28 
respectively. Suppose *beta0* changes to values [100, 103, 107, 200] at steps [2, 15, 17, 
20] and *mu* changes to values [0.1, 0.3, 0.8] at steps [ 6, 15, 20]. Then,
 - *breaks* = [2, 6, 15, 17, 20]
 - *count_of_tv_at_breaks* = [1, 1, 2, 1, 2]
 - *tv_spi* = [15, 28, 15, 28, 15, 15, 28]
 - *tv_val* = [100, 0.1, 103, 0.3, 107, 200, 0.8]

## Update rate matrix
Based on the above data structure, the algorithm of *update_ratemat* will have two phases: 
 1. Update vector *sp* to make sure every time-varying parameter get updated by each step 
 of simulation. 
 2. Loop through *updateidx* to update the rate matrix.

Below is a C-like pseudo code:

     int n = 300; // number of simulations steps int nextBreak = 
    0; int start = 0; for (i=0; i<n; i++) {
        // 1 update sp
        if (nextBreak<breaks.size() && i==breaks[nextBreak]) { for (j=start; 
            j<start+count_of_tv_at_break[nextBreak]; j++)
                sp[tv_spi[j]] = tv_val[j]; start += count_of_tv_at_break[nextBreak]; 
            nextBreak++;
        }
        // 2 Update those elements that need to be updated
        for (j=0; j<updateidx.size(); j++) { Update element [from[updateidx[j]], 
            to[updateidx[j]] by using formula defined by count[updateidx[j]] plus 
            corresponding spi and modifier, and data sp.
        }
    }

## Future considerations
### Considerations about spec 0.0.4
We can use the above structure to cover exogenous time-varying parameters by pre-computing 
its values at all the time stpes and treating it the same way as we do for the piece-wise 
parameters. However, we can design a better way to handle exogenous time-varying 
parameters.

### Pack up parameters
In previous implementations of TMB/C++, we pass all the elements in the data structure 
separately. A better way is to pack all these elements (such as *sp*, *from*, *to*, 
*count*, *spi*, *modifier*, *breaks*, etc) into a TMB's 
[DATA_STRUCT](https://kaskr.github.io/adcomp/group__macros.html#gaf9885566da0d248c1a4b4d7a0eeafcd2). 
Perhaps, we can pack them into hierarchical struct, e.g.,
 - Packing *from*, *to* and *count* into *non_elements";
 - Packing *spi* and *modifier* into *operands*;
 - Packing *non_elements" and *operands* into *formulas* or *ratemat*

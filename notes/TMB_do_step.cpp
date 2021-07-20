#include <iostream>
#include <TMB.hpp>

//#define TMB_VERBOSE
//using namespace std;	// causes trouble with definitions in TMB.hpp
static int count = 0;

// simplest/stripped-down do_step
// no hazard, no stochasticity, no 'exponential step'
template<class Type>
Type objective_function<Type>::operator() ()
{
    std::cout << "================== Begin of TMB_do_step() =================" << std::endl;

    std::cout << "static variable = " << count++ << std::endl;

    // Joint negative log-likelihood (stub)
    Type jnll=0;

    DATA_VECTOR(state);
    DATA_SPARSE_MATRIX(ratemat);
    DATA_IVECTOR(inf_ind);       // indices of infectious classes
    DATA_INTEGER(transm_ind);    // index of contact rate
    DATA_IVECTOR(transm_wt_ind); // index of contact rate modifiers
    DATA_IVECTOR(foi_ind);       // indices (i,j) of foi element of ratemat

    // 1-based indexing in R to 0-based indexing in C++
    // DATA_IVECTOR is defined as
    //   template <class Type>
    //   struct vector : Array<Type,Dynamic,1>
    // where Type is int
    // The sentence below doesn't work.
    //std::transform(inf_ind.begin(), inf_ind.end(), inf_ind.begin(), [](int x) -> int {return (x-1);});

#ifdef TMB_VERBOSE
    std::cout << "state = " << state << std::endl;
    std::cout << "ratemat = " << ratemat << std::endl;
    std::cout << "inf_ind = "<< inf_ind << std::endl;
    std::cout << "transm_ind = " << transm_ind << std::endl;
    std::cout << "transm_wt_ind = " << transm_wt_ind << std::endl;
    std::cout << "foi_ind = " << foi_ind << std::endl;
#endif
    // eventually we will want gradients with respect to these parameters (although by that time this will be a function within
    // a larger objective function)
    // PARAMETER_* denotes an object whose gradients are propagated
    PARAMETER_VECTOR(params);

#ifdef TMB_VERBOSE
    std::cout << "params = " << params << std::endl;
#endif

    int ns = state.size();
   
    // update force of infection (foi) in rate matrix

    Type foi = 0;
    for (int i = 0; i < inf_ind.size(); i++) {
        // note TMB vectors use () not [] for indexing (but still zero-indexed)
        foi += state(inf_ind(i))*params(transm_wt_ind(i));
    }
    foi *= params(transm_ind);

    // need special accessors for Eigen (sparse) rate matrix operators:
    // https://stackoverflow.com/questions/42376127/how-to-access-a-specific-row-col-index-in-an-c-eigen-sparse-matrix
    // don't know if this operation will work correctly for AD propagation? (i.e., we're updating a DATA_* type object ...)
    ratemat.coeffRef(foi_ind(0), foi_ind(1)) = foi;
   
    // transform ratemat to flows
    // for now, we can transform; later we might need to copy/update independently, if we are updating ratemat across steps
    // should take better advantage of sparsity!
    // does Eigen and/or TMB have a built-in column-multiply operator?
    // might eventually want to have persistent/pre-allocated ratemat *and* flowmat objects
    for (int i=0; i < state.size(); i++) {
        for (int j=0; j < state.size(); j++) {
            if (ratemat.coeff(i,j) != 0.0) {
	        ratemat.coeffRef(i,j) = ratemat.coeff(i,j) * state(j);
	    }
        }
    }
   
    // calculate outflows & inflows
    vector<Type> outflow(ns);
    for (int i=0; i < ns; i++) {
        outflow(i) = ratemat.row(i).sum();
    }

    vector<Type> inflow(ns);
    for (int i=0; i < ns; i++) {
        inflow(i) = ratemat.col(i).sum();
    }

    // update state
    state = state - outflow + inflow;

    // 'DATA_*' objects are _not_ automatically mutable
    // may use DATA_UPDATE() to propagate changes back to R

    // makes state available 
    REPORT(state);

    std::cout << "================== End of TMB_do_step() =================" << std::endl;
    std::cout << std::flush;

    jnll = 777;
    return jnll;
}
   



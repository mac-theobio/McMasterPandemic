#include <iostream>
#include <sys/time.h>
#include <TMB.hpp>

//#define TMB_VERBOSE

template<class Type>
Type objective_function<Type>::operator() ()
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    std::cout << "*** Begin of TMB_do_step() at " << tv.tv_sec*1000000+tv.tv_usec << std::endl;

    // Joint negative log-likelihood (stub)
    Type jnll=789;

    DATA_VECTOR(state);
    DATA_SPARSE_MATRIX(ratemat);
    DATA_INTEGER(dt);			// integer->integer
    DATA_INTEGER(do_hazard);		// boolean->integer
    DATA_INTEGER(stoch_proc);		// boolean->integer
    DATA_INTEGER(do_exponential);	// boolean->integer
    DATA_STRING(testwt_scale);		// string->string

#ifdef TMB_VERBOSE
    std::cout << "state = " << state << std::endl;
    std::cout << "ratemat = " << ratemat << std::endl;
    std::cout << "dt = "<< dt << std::endl;
    std::cout << "do_hazard = " << do_hazard << std::endl;
    std::cout << "stoch_proc = " << stoch_proc << std::endl;
    std::cout << "do_exponential = " << do_exponential << std::endl;
    std::cout << "testwt_scale = "<< testwt_scale << std::endl;
#endif
    // eventually we will want gradients with respect to these parameters (although by that time this will be a function within
    // a larger objective function)
    // PARAMETER_* denotes an object whose gradients are propagated
    PARAMETER_VECTOR(params);

#ifdef TMB_VERBOSE
    std::cout << "params = " << params << std::endl;
#endif

    gettimeofday(&tv, NULL);
    std::cout << "*** End of TMB_do_step() at   " << tv.tv_sec*1000000+tv.tv_usec << std::endl;

    std::cout << std::flush;

    return jnll;
}
   



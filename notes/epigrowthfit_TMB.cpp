#include <TMB.hpp>
#include <iostream>
#include <string>

// eventual goal: flexible epigrowthfit-style fitting
// possibility of including multiple time series with random effects parameters
// compare cumulative vs incidence fitting?
// regularization? priors?

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}



// https://kaskr.github.io/adcomp/_book/Tutorial.html
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(t);      // times (includes step [-1])
  DATA_VECTOR(x);      // initial densities
  
  // Parameters  
  PARAMETER(log_thalf);    // log half-max time
  PARAMETER(log_K);        // log carrying capacity
  PARAMETER(log_r);        // log growth rate
  PARAMETER(log_p);        // log Richards shape parameter
  PARAMETER(log_nb_disp);  // log NB dispersion parameter (FIXME: reparameterize?)

  vector<Type> cumcurve(t.size()); // includes starting time
  vector<Type> inccurve(x.size());
  
  // Objective funcction
  Type jnll = 0;

  // inverse-link functions
  Type thalf = exp(log_thalf);
  Type K = exp(log_K);
  Type r = exp(log_r);
  Type p = exp(log_p);
  Type nb_disp = exp(log_nb_disp); 
  
  // parameter est

  for (int i = 0; i < t.size(); i++) {
    cumcurve(i) = K/(1+exp(-r*(t(i)-thalf)));
  }
  for (int i = 0; i < x.size(); i++) {
    inccurve(i) = cumcurve(i+1)-cumcurve(i);
    if (!isNA(x(i))) {
      // std::cout << i << " " << x(i) << " " << inccurve(i)  << std::endl;      
      jnll -= dnbinom2(x(i),inccurve(i),nb_disp,1);
    }
  }
  
  REPORT(inccurve);
  return jnll;
}

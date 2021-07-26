#include <TMB.hpp>

// simple SIR model: trajectory-matching model
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Joint negative log-likelihood (stub)
  Type jnll=0;

   DATA_VECTOR(log_obs_incidence);
   DATA_SCALAR(N);
	
   PARAMETER(log_beta);
   PARAMETER(log_gamma);
   PARAMETER(log_I0);
   PARAMETER(log_nbdisp);

   Type S, log_I;

   int nobs = log_obs_incidence.size();
   vector<Type> log_incidence(nobs);

   Type I0 = exp(log_I0);
   S = N - I0;
   log_I = log_I0;
   
   for (int t = 0; t < nobs ; t++) {
      log_incidence(t) = log_beta + log_I + log(S);
      S -= exp(log_incidence(t));
      log_I += exp(log_beta) - exp(log_gamma);
      // dnbinom_robust parameterized by log(mu), log(var-mu)
      jnll -= dnbinom_robust(log_obs_incidence(t), log_incidence(t),
			     2. * log_incidence(t) - log_nbdisp, true);
      SIMULATE {
	      Type s1 = exp(log_incidence(t));      
	      Type s2 = s1 * (Type(1) + s1 / exp(log_nbdisp));
	      log_obs_incidence(t) = rnbinom2(s1, s2);
	      
      }
   }
   
   return jnll;
}
   



#include <TMB.hpp>

// simple SIR model: trajectory-matching model
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Joint negative log-likelihood
  Type jnll=0;

  DATA_VECTOR(obs_incidence);      // observed incidence (new cases) over time
  DATA_SCALAR(N);                  // pop size

  // chose to treat parameters as scalars rather than elements of a parameter vector
  // do_step.cpp indexes the parameter vector in order to allow for (lots) of flexibility
  //   in the types of models that can be handled
  // fit positive parameters on log scale; this is handled via 'link functions' in MacPan
  PARAMETER(log_beta);    // transmission parameter ('contact rate')
  PARAMETER(log_gamma);   // recovery rate
  //PARAMETER(log_I0);      // initial number infected (proportion?)
  PARAMETER(log_nbdisp);  // negative binomial dispersion parameter

  Type S, log_I;  // state variables (scalars)

  int nobs = obs_incidence.size();  // number of observations
  Type log_I0 = -3; // FIXME: remove hardcoding
  vector<Type> log_incidence(nobs); // model-based/theoretical incidence

   Type I0 = exp(log_I0);
   S = N - I0;      // the susceptible population is the total (N) minus
                    // the initial number infected (I0)
   log_I = log_I0;  // evaluate I dynamics on the log scale; I is the most
                    // likely value to get in trouble with dipping negative
                    // this is *not* done in MacPan because (?) it doesn't fit
                    // nicely into the overall 'flows' framework

   for (int t = 0; t < nobs ; t++) {
      log_incidence(t) = log_beta + log_I + log(S);
      // no delta-t here yet ... time steps are dt=1
      S -= exp(log_incidence(t));
      log_I += exp(log_beta) - exp(log_gamma);  // d(log(I))/dt = (dI/dt)/I = beta-gamma

      // this part is a little bit mind-bending
      // TMB has two parameterizations of the negative binomial likelihood
      //   one is the (prob, 'size') from math stats (i.e. the R default)
      //   the other is 'robust'; it takes log(mean) and log(var/mu)
      //   ('robust' because if we're doing computations of the mean on the log scale
      //    we don't have to back-transform and lose precision in order to calculate the
      //    log-likelihood)
      // R's alternative is (mu, size) where var/mu = 1+mu/size
      // dnbinom_robust parameterized by log(mu), log(var/mu)
      // I used this one because I cheated and copied code from glmmTMB
      jnll -= dnbinom_robust(obs_incidence(t),   // observed value
			     log_incidence(t),   // log(predicted mean)
			     2. * log_incidence(t) - log_nbdisp, // stolen from glmmTMB w/o checking
			     true                // we want the log-likelihood
			     );
      SIMULATE {
	      // simulation code: adapted from glmmTMB
	      Type s1 = exp(log_incidence(t));
	      Type s2 = s1 * (Type(1) + s1 / exp(log_nbdisp));
	      obs_incidence(t) = rnbinom2(s1, s2);
	      REPORT(obs_incidence);
	      REPORT(log_incidence);
      }
   }

   return jnll;
}

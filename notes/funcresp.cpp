#include <TMB.hpp> // for print statements
#include <iostream>
#include <string>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Step 1: Code up a plain C version (give a 1 number (double), and it will output 1 number (double))
double W(double x) { // defining types for argument in TMB EAH: (takes and returns)
  double reltol = 1.490116e-08; // define tolerance for iterations, smallest x s.t. 1+x != 1 EAH: double=numeric, tol in R: sqrt(.Machine$double.eps), smallest x s.t. 1+x != 1
  double logx = log(x); 
  double y = (logx > 0 ? logx : 0); // EAH: C++ -  if first is true, return 2nd, otherwise do 3rd
  int niter = 100; // maximum number of iterations
  int i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < reltol) break; // solve y=e^y, with newton's method : could substitute faster C++ solver here
    double fdivg = (y - exp(logx - y)) / (1 + y);
    y = y - fdivg;
  }
  if (i == niter) Rf_warning("W: failed convergence"); // return warning if you hit the # of iterations
  return y; //return value
}

// Step 2: Code up reverse mode derivatives:
// x = W(x) * exp( W(x) )
// W'(x) = 1 / ( exp(W(x)) * (1 + W(x)) ) (derivation on Wikipedia)

TMB_ATOMIC_VECTOR_FUNCTION( // tells TMB what the derivative of the new function is
  // ATOMIC_NAME
  W
  ,
  // OUTPUT_DIM 
  1,//(EAH:1=vector , separate)
  // ATOMIC_DOUBLE
  double x = tx[0]; // input (defining the input for the wrapper function)
ty[0] = W(x);     // output
,
// ATOMIC_REVERSE
Type W = ty[0];   // Function value (defining the derivative, need dy/dx in terms of y (reverse mode)) type so it's flexible for shape
px[0] = 1. / (exp(W) * (1. + W)) * py[0]; //using chain rule for the derivative (line 21), whole thing is dW/dy * dy/dx, py[0] is dy/dx
) // end of TMB_ATOMIC_VECTOR_FUNCTION
  
  // Step 3: Make version with scalar arguments: EAH:  More wrapping - making types more general than a double
  // (input number of any type, output number of the same type as the input)
  template<class Type>
  Type W(Type x){
    CppAD::vector<Type> tx(1);//EAH: what does tx(1): length (i.e. a scalar)
    tx[0] = x;
    return W(tx)[0];
  }

// see glmmTMB code for enum handling between C++ and R, in particular
//  .../R/enum.r (see Makefile for automated generation of enum.R)
// source("enum.R")
enum size_model_type {
  constant  = 0,
  ricker    = 1,
  powricker = 2
};

enum consume_model_type {
  holling = 0,
  rogers = 1
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_IVECTOR(block);   // experimental blocks
  DATA_VECTOR(killed);  // number of prey killed              
  DATA_VECTOR(N0);      // initial densities
  DATA_VECTOR(Sizes);   // size classes
  DATA_INTEGER(ndata);  // number of data rows
  DATA_INTEGER(nblock); // number of blocks (random effects)
  DATA_INTEGER(tval);  // length of experiment (number of timesteps)
  DATA_INTEGER(size_model_flag);    // integer code for size model
  DATA_INTEGER(consume_model_flag); // integer code for consumption model
  
  // Parameters
  
  PARAMETER_VECTOR(consume_parms); // parameter vector
  PARAMETER_VECTOR(cRE);// random effects
  PARAMETER(log_SD_cRE);// random effect variance
  
  Type g; // gamma for power-Ricker
  Type d; // size at max attack rate
  Type h; // handling time
  Type c; // max attack rate
  Type a; // size-dependent attack rates
  Type jnll = 0;
  Type risk;
  
  int p=0;  // parameter counter
  if(size_model_flag==constant || size_model_flag==ricker || size_model_flag==powricker) {
    h = exp(consume_parms(p++)); // handling time
    c = exp(consume_parms(p++)); // max attack rate
    // logc = consume_parms(p++);   // then exponentiate later
  }
  
  if(size_model_flag==ricker || size_model_flag==powricker) {
    d = exp(consume_parms(p++)); // size at max attack rate
  }
  
  if(size_model_flag==powricker){
    g = exp(consume_parms(p++)); // shape of power ricker curve
  }
  
  // std::cout << "size_model_flag " << size_model_flag << " " << constant << std::endl;
  // parameter est
  for (int i=0; i<ndata; i++) {
    switch (size_model_flag) {
    case constant:
      a = c+cRE(block(i));  //adding random effect on linear scale
      // or, a = exp(logc +cRE(block(i)));
      break;
    case ricker:
      // Probability of data conditional on fixed and random effect values
      a = (c+cRE(block(i)))*(Sizes(i)/d)*exp(1-Sizes(i)/d); // Ricker size-dependent attack rate
      // std::cout << "size loop " << a << " " << std::endl;
      // or
      // loga = (logc + cRE(block(i))) + log(Sizes(i)) -logd +(1-Sizes(i)/d)
    case powricker:
      a = (c+cRE(block(i)))*pow(Sizes(i)/d,g)*(exp(1-Sizes(i)/d));
      // loga = (logc + cRE(block(i))) + g*(log(Sizes(i)) -logd) +(1-Sizes(i)/d)
      break;
    default:
      error("unknown size model!");
    } // End size switch
    // std::cout << "consume_model_flag " << consume_model_flag << " " << holling << std::endl;


    // compute predation risk
    switch (consume_model_flag) {
    case holling:
      // std::cout << "risk calc " << i << " " << a << " " << h << " " << N0(i) << std::endl;
      risk = a/(1+a*h*N0(i));
      // std::cout << "risk loop " << risk << std::endl;
      break;
    case rogers:
      // Time set to 1 for now
      risk = N0(i) - W(a*h*N0(i)*exp(-a*(tval-h*N0(i))))/(a*h);
    // risk = N0(i) - W(a*h*tval); // Rogers eq (copied from FR Tutorial), now we need time... setting as 1 right now
      risk /= N0(i);
      // FIXME: could express risk directly as
      //   1- W(...)/(a*h*N0(i));
      break;
    default:
      error("unknown consumption model!");
      //std::cout << risk << std::endl;
    } // end consumption model switch

    // FIXME: allow beta-binomial? Could steal it from glmmTMB.cpp
    // std::cout << "loglik calc" << i << " " << killed(i) << " " << N0(i) << " " << risk << std::endl;

    Type nlli = -dbinom(killed(i),N0(i),risk,true);
    std::cout << "i " << c << " " << cRE(block(i)) << " " << Sizes(i) << " " << d << " " << a << " " << h << " " << risk << " " << killed(i) << " " << N0(i) << " " << nlli << std::endl;
    if ( !isNA(nlli)) {
      jnll += nlli;
      std::cout << "OK" << std::endl;
    } else {
      jnll += 99999;
      std::cout << "bad" << std::endl;
    }
  } // end data loop
  // Probability of random coefficients
  for( int i=0; i<nblock; i++){
    jnll += -dnorm(cRE(i), Type(0.0), exp(log_SD_cRE), true);
  }

  // Reporting
  ADREPORT(cRE);

  //std::cout << "finished loop" << std::endl;
  //jnll += -sum(dnorm(logc, Type(0.0), sigmaC, true)); // RE term for block effect
  //ADREPORT(logSigma);
  
  return jnll;
}

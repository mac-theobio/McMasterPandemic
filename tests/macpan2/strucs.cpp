// Simple linear regression.
#include <Eigen/Eigen>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h> 	// isnan() is defined
#include <sys/time.h>
#include <TMB.hpp>
#include <cppad/local/cond_exp.hpp>


template<class Type>
struct ll {
  // below is a vector of vectors that passed from R
  vector<vector<Type>> vectors;
  
  ll(SEXP ii){ // Constructor
    // Get elements by their names
    //Y = asVector<Type>(getListElement(ii,"Y"));
    //x = asVector<Type>(getListElement(ii,"x"));

    // Get elements by their indices
    int n = length(ii);
    vector<vector<Type>> vs(n);
    vectors = vs;

    for (int i = 0; i < n; i++) 
      vectors[i] = asVector<Type>(VECTOR_ELT(ii, i));
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{

  //SEXP MatrixList = getListElement(TMB_OBJECTIVE_PTR -> data, "MatrixList")

  //DATA_VECTOR(Y);
  //DATA_VECTOR(x);
  DATA_STRUCT(Yx, ll);
  for (int i=0; i<Yx.vectors.size(); i++) {
    std::cout << "i = " << i << std::endl;
    std::cout << "vector = " << Yx.vectors[i] << std::endl;
  }

  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));
  Type nll = -sum(dnorm(Yx.vectors[0], a+b*Yx.vectors[1], exp(logSigma), true));
  return nll;
}

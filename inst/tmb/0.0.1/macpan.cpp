// update_ratemat2 is the function that applies our own "parser" to the formulas
// to calculate the non-zero elements in rate matrix.

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/time.h>
#include <TMB.hpp>

///////////////////////////////////////////////////////////////////////////////
// Helper function update_ratemat2
template<class Type>
Eigen::SparseMatrix<Type> update_ratemat2(
    const Eigen::SparseMatrix<Type>& ratemat,
    const vector<Type>& sp,
    const vector<int>& from,
    const vector<int>& to,
    const vector<int>& count,
    const vector<int>& spi,
    const vector<int>& modifier)
{
  Eigen::SparseMatrix<Type> result(ratemat.rows(), ratemat.cols());

  //std::cout << "=========================================" << std::endl;
  //std::cout << ratemat << std::endl;

  int start = 0;
  int n = count.size();
  for (int i=0; i<n; i++) {
    int row = from[i] - 1;
    int col = to[i] - 1;
    result.coeffRef(row,col) = 0.0;
    Type prod = 1.0;
    for (int j=start; j<start+count[i]; j++) {
      Type x = sp[spi[j]-1];
      if (modifier[j] & 0b100) {
        result.coeffRef(row,col) += prod;
        prod = 1;
      }
      if (modifier[j] & 0b001)
        x = 1-x;
      else if (modifier[j] & 0b010)
        x = 1/x;
      prod *= x;
    }
    result.coeffRef(row,col) += prod;
    start += count[i];
  }

  //std::cout << "-------------------------------------" << std::endl;
  //std::cout << result << std::endl;
  return result;
}


///////////////////////////////////////////////////////////////////////////////
// Helper function col_multiply
template<class Type>
Eigen::SparseMatrix<Type> col_multiply(
    const Eigen::SparseMatrix<Type>& mat,
    const vector<Type>& vec)
{
  Eigen::SparseMatrix<Type> newmat = mat;
  Type* valPtr = newmat.valuePtr();
  int* rowIndexPtr = newmat.innerIndexPtr();
  int* outPtr = newmat.outerIndexPtr();

  for (int j= 0; j<newmat.outerSize(); j++)
    for (int i= outPtr[j]; i<outPtr[j+1]; i++) {
      int row = rowIndexPtr[i];
      valPtr[i] = valPtr[i]*vec(row);
    }

    return newmat;
}

///////////////////////////////////////////////////////////////////////////////
// Helper function rowSums
template<class Type>
vector<Type> rowSums(
    const Eigen::SparseMatrix<Type>& mat)
{
  vector<Type> result(mat.rows());

  for (int i=0; i< mat.rows(); i++)
    result(i) = mat.row(i).sum();

  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Helper function colSums
template<class Type>
vector<Type> colSums(
    const Eigen::SparseMatrix<Type>& mat)
{
  vector<Type> result(mat.cols());

  for (int i=0; i< mat.cols(); i++)
    result(i) = mat.col(i).sum();

  return result;
}

///////////////////////////////////////////////////////////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{

  // Joint negative log-likelihood (stub)
  Type jnll= 0;

  // Get data and parameters from R
  DATA_VECTOR(state);
  DATA_SPARSE_MATRIX(ratemat);
  DATA_IVECTOR(from);
  DATA_IVECTOR(to);
  DATA_IVECTOR(count);
  DATA_IVECTOR(spi);
  DATA_IVECTOR(modifier);
  //DATA_INTEGER(dt);			// integer->integer
  //DATA_INTEGER(do_hazard);		// boolean->integer
  //DATA_INTEGER(stoch_proc);		// boolean->integer
  //DATA_INTEGER(do_exponential);	// boolean->integer
  //DATA_STRING(testwt_scale);		// string->string

  PARAMETER_VECTOR(params);

  // Concatenate state and params
  vector<Type> sp(state.size()+params.size());
  sp << state, params;
  std::cout << "using 0.0.1" << std::endl;

  // We've got everything we need, lets do the job ...
  ratemat = update_ratemat2(ratemat, sp, from, to, count, spi, modifier);

  REPORT(ratemat);

  return jnll;
}

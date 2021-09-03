// update_ratemat2 is the function that applies our own "parser" to the formulas
// to calculate the non-zero elements in rate matrix.

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/time.h>
#include <TMB.hpp>

///////////////////////////////////////////////////////////////////////////////
// This function doesn't remove the columns. Instead, it put zeros for these
// columns.
template<class Type>
void remove_cols(
    Eigen::SparseMatrix<Type>& mat,
    const vector<int>& indices_to_remove)
{
  //std::cout << "matrix before= " << mat << std::endl;

  Type* valPtr = mat.valuePtr();
  int* rowIndexPtr = mat.innerIndexPtr();
  int* outPtr = mat.outerIndexPtr();

  for (int j= 0; j<indices_to_remove.size(); j++) {
    int jj = indices_to_remove[j] - 1;
    //std::cout << "jj= " << jj << std::endl;
    for (int i= outPtr[jj]; i<outPtr[jj+1]; i++) {
      //std::cout << valPtr[i] << " will be zero-ed out ****" << std::endl;
      valPtr[i] = 0.0;
    }
  }
  //std::cout << "matrix after= " << mat << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Helper function make_ratemat
template<class Type>
Eigen::SparseMatrix<Type> make_ratemat(
    int size,
    const vector<Type>& sp,
    const vector<int>& from,
    const vector<int>& to,
    const vector<int>& count,
    const vector<int>& spi,
    const vector<int>& modifier)
{
  Eigen::SparseMatrix<Type> result(size, size);
  //result.insert(2,3) = 23;
  //result.insert(4,13) = 413;
  //result.insert(9,13) = 913;
  //result.makeCompressed();
  //remove_cols(result, from);

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
  
  result.makeCompressed();

  //std::cout << "-------------------------------------" << std::endl;
  //std::cout << result << std::endl;
  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Helper function update_ratemat
template<class Type>
Eigen::SparseMatrix<Type> update_ratemat(
    const Eigen::SparseMatrix<Type>& ratemat,
    const vector<Type>& sp,
    const vector<int>& from,
    const vector<int>& to,
    const vector<int>& count,
    const vector<int>& spi,
    const vector<int>& modifier)
{
  Eigen::SparseMatrix<Type> result = ratemat; //(ratemat.rows(), ratemat.cols());

  //std::cout << "=========================================" << std::endl;
  //std::cout << ratemat << std::endl;

  int start = 0;
  int n = count.size();
  //std::cout << "n = " << n << std::endl;

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
  std::cout << "===================== objective_function ====================" << std::endl;

  // Joint negative log-likelihood (stub)
  Type jnll= 0;

  // Get data and parameters from R
  DATA_VECTOR(state);
  //DATA_SPARSE_MATRIX(ratemat);
  DATA_IVECTOR(from);
  DATA_IVECTOR(to);
  DATA_IVECTOR(count);
  DATA_IVECTOR(spi);
  DATA_IVECTOR(modifier);
  DATA_IVECTOR(update_from);
  DATA_IVECTOR(update_to);
  DATA_IVECTOR(update_count);
  DATA_IVECTOR(update_spi);
  DATA_IVECTOR(update_modifier);
  DATA_IVECTOR(par_accum_indices);

  //std::cout << "update_from = " << update_from << std::endl;
  //std::cout << "update_count = " << update_count << std::endl;
  //std::cout << "update_spi = " << update_spi << std::endl;
  //std::cout << "par_accum_indices = " << par_accum_indices << std::endl;

  PARAMETER_VECTOR(params);

  // Concatenate state and params
  vector<Type> sp(state.size()+params.size());
  sp << state, params;
  //std::cout << "sp = " << sp << std::endl;

  // We've got everything we need, lets do the job ...
  Eigen::SparseMatrix<Type> ratemat = make_ratemat(state.size(), sp, from, to, count, spi, modifier);

  int numIterations = 3; //10000;
  int stateSize = state.size();
  vector<Type> concatenated_state_vector(numIterations*stateSize);
  //concatenated_state_vector.block(0, 0, stateSize, 1) = state;
  for (int i=0; i<numIterations; i++) {
    ratemat = update_ratemat(ratemat, sp, update_from, update_to, update_count, update_spi, update_modifier);
    Eigen::SparseMatrix<Type> flows = col_multiply(ratemat, state); // ignore * dt
    vector<Type> inflow = colSums(flows);
    remove_cols(flows, par_accum_indices);
    vector<Type> outflow = rowSums(flows); // remove some columns before doing so
    state = state - outflow + inflow;
    sp << state, params;
    concatenated_state_vector.block(i*stateSize, 0, stateSize, 1) = state;
  }

  //std::cout << "concatenated_state_vector = " << concatenated_state_vector << std::endl;

  REPORT(ratemat);
  REPORT(concatenated_state_vector);

  return jnll;
}

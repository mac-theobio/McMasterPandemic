// This version implements spec 0.1.1 in https://canmod.net/misc/flex_specs

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/time.h>
#include <TMB.hpp>


///////////////////////////////////////////////////////////////////////////////
// Helper function rowSums
template<class Type>
vector<Type> rowSums(
    const Eigen::SparseMatrix<Type>& mat)
{
  vector<Type> result(mat.rows());

  for (int i=0; i< mat.rows(); i++)
    // here we need to zero out flows
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
// This function doesn't remove the columns. Instead, it put zeros for these
// columns.
template<class Type>
void remove_cols(
    Eigen::SparseMatrix<Type>& mat,
    const vector<int>& indices_to_remove)
{
  Type* valPtr = mat.valuePtr();
  int* outPtr = mat.outerIndexPtr();

  // loop over columns to zero-out
  for (int j= 0; j<indices_to_remove.size(); j++) {
    int jj = indices_to_remove[j] - 1;

    // loop over all rows in the current column
    for (int i= outPtr[jj]; i<outPtr[jj+1]; i++) {
      valPtr[i] = 0.0;
    }
  }
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
        if (x > 1e-5) {
          x = 1/x;
        }
      prod *= x;
    }
    result.coeffRef(row,col) += prod;
    start += count[i];
  }

  result.makeCompressed();

  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Helper function update_ratemat
template<class Type>
void update_ratemat(
    Eigen::SparseMatrix<Type>* ratemat, // to be updated
    const vector<Type>& sp,
    const vector<int>& from,
    const vector<int>& to,
    const vector<int>& count_integral,
    const vector<int>& spi,
    const vector<int>& modifier,
    const vector<int>& updateidx)
{
  for (int i=0; i<updateidx.size(); i++) {
    int idx = updateidx[i] - 1;
    int row = from[idx] - 1;
    int col = to[idx] - 1;

    ratemat->coeffRef(row,col) = 0.0;
    Type prod = 1.0;
    for (int j=count_integral[idx]; j<count_integral[idx+1]; j++) {
      Type x = sp[spi[j]-1];
      if (modifier[j] & 0b100) {
        ratemat->coeffRef(row,col) += prod;
        prod = 1;
      }
      if (modifier[j] & 0b001)
        x = 1-x;
      else if (modifier[j] & 0b010)
        if (x > 1e-5) {
          x = 1/x;
        }
      prod *= x;
    }
    ratemat->coeffRef(row,col) += prod;
  }
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
// Helper function calc_flowmat
template<class Type>
Eigen::SparseMatrix<Type> calc_flowmat(
    const Eigen::SparseMatrix<Type>& mat,
    const vector<Type>& vec,
    int do_hazard)
{
  if (!do_hazard)
    return col_multiply(mat, vec);
  else {
    vector<Type> r = rowSums(mat);
    vector<Type> rho = exp(-r);

    vector<Type> s_tilde(vec.size());
    for (int i=0; i<vec.size(); i++)
      if (r[i]==0) // shall it be something like "abs(r[i])<0.00001" ?
        s_tilde[i] = 0.0;
      else
        s_tilde[i] = vec[i]/r[i];

    vector<Type> v = s_tilde*(1.0-rho);

    // Maybe a more efficient way is using a modified version of col_multiply()
    //Eigen::SparseMatrix<Type> flowmat = col_multiply(mat, v);
    //flowmat.diagonal() = 0.0; // doesn't work
    //return flowmat;

    Eigen::SparseMatrix<Type> newmat = mat;
    Type* valPtr = newmat.valuePtr();
    int* rowIndexPtr = newmat.innerIndexPtr();
    int* outPtr = newmat.outerIndexPtr();

    for (int j= 0; j<newmat.outerSize(); j++)
      for (int i= outPtr[j]; i<outPtr[j+1]; i++) {
        int row = rowIndexPtr[i];
        if (row==j)
           valPtr[i] = 0.0;
        else
          valPtr[i] = valPtr[i]*v[row];
      }

    return newmat;
  }
}

///////////////////////////////////////////////////////////////////////////////
// Helper function update_sum_in_sp
template<class Type>
void update_sum_in_sp(
    vector<Type>& sp,
    const vector<int>& sumidx,
    const vector<int>& sumcount,
    const vector<int>& summandidx)
{
  int start = 0;
  for (int i=0; i<sumidx.size(); i++) {
    int idx = sumidx[i] - 1;
    sp[idx] = 0.0;

    for (int j=start; j<start+sumcount[i]; j++) {
      sp[idx] += sp[summandidx[j]-1];
    }
    start += sumcount[i];
    //std::cout << "idx = " << idx << std::endl;
    //std::cout << "sp = " << sp[idx] << std::endl;
  }
}

template<class Type>
vector<Type> do_step(
    vector<Type> state,
    Eigen::SparseMatrix<Type> ratemat,
    vector<int> par_accum_indices,
    int do_hazard
)
{
  // Calculate flow matrix
  Eigen::SparseMatrix<Type> flows = calc_flowmat(ratemat, state, do_hazard);
  vector<Type> inflow = colSums(flows);
  remove_cols(flows, par_accum_indices);
  vector<Type> outflow = rowSums(flows); // remove some columns before doing so
  state = state - outflow + inflow;
  return state;
}

template <class Type>
struct update_state_functor{

  Eigen::SparseMatrix<Type> ratemat_;
  vector<int> par_accum_indices_;
  int do_hazard_;

  // constructor
  update_state_functor(
    Eigen::SparseMatrix<Type> ratemat,
    vector<int> par_accum_indices,
    int do_hazard) {
      std::cout << "here in the constructor...";
      ratemat_ = ratemat;
      par_accum_indices_ = par_accum_indices;
      do_hazard_ = do_hazard;

  }

  template <typename T>
  vector<T> operator()(vector<T> state_) {
    // 1 transform state from vector<T> to vector<Type>
    int n = state_.size();
    vector<Type> st(n);
    for(int i=0; i<n; i++)
       st[i] = CppAD::Value((AD<Type>)state_[i]);

    // 2 do all the calculations in Type
    vector<Type> updated_state = do_step(st, ratemat_, par_accum_indices_, do_hazard_);

    // 3 transform final result from vector<Type> back to vector<T>
    CppAD::vector<T> xx = CppAD::vector<T>(updated_state);
    std::cout << "here in the functor..." << updated_state.coeff(0) << "..." << xx[0];
    //return updated_state;
    return (xx);
  }

};

///////////////////////////////////////////////////////////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Joint negative log-likelihood (stub)
  //Type jnll= 0;

  // Get data and parameters from R
  DATA_VECTOR(state);
  //DATA_SPARSE_MATRIX(ratemat);
  DATA_IVECTOR(from);
  DATA_IVECTOR(to);
  DATA_IVECTOR(count);
  DATA_IVECTOR(spi);
  DATA_IVECTOR(modifier);
  DATA_IVECTOR(updateidx);
  DATA_IVECTOR(breaks);
  DATA_IVECTOR(count_of_tv_at_breaks);
  DATA_IVECTOR(tv_spi);
  //DATA_VECTOR(tv_val);
  DATA_VECTOR(tv_mult);
  DATA_IVECTOR(tv_orig);
  DATA_IVECTOR(par_accum_indices);

  DATA_IVECTOR(sumidx);
  DATA_IVECTOR(sumcount);
  DATA_IVECTOR(summandidx);

  DATA_INTEGER(numIterations);
  DATA_INTEGER(do_hazard);

  PARAMETER_VECTOR(params);

  std::cout << "here in the objective function...";
  // state = make_state(params);

  //std::cout << "breaks = " << breaks << std::endl;
  //std::cout << "count_of_tv_at_breaks = " << count_of_tv_at_breaks << std::endl;
  //std::cout << "tv_spi = " << tv_spi << std::endl;
  //std::cout << "tv_val = " << tv_val << std::endl;

  //std::cout << "tv_mult = " << tv_mult << std::endl;
  //std::cout << "tv_orig = " << tv_orig << std::endl;

  // Concatenate state and params
  vector<Type> sp(state.size()+params.size()+sumidx.size());
  vector<Type> place_holder(sumidx.size());
  sp << state, params, place_holder;

  update_sum_in_sp(sp, sumidx, sumcount, summandidx);

  // make a copy of sp
  vector<Type> sp_orig(sp);

  //std::cout << "sp = " << sp << std::endl;
  //std::cout << "sp_orig = " << sp_orig << std::endl;

  // Calculate integral of count
  vector<int> count_integral(count.size()+1);
  count_integral[0] = 0;
  for (int i=0; i<count.size(); i++)
    count_integral[i+1] = count_integral[i] + count[i];

  // We've got everything we need, lets do the job ...
  Eigen::SparseMatrix<Type> ratemat = make_ratemat(state.size(), sp, from, to, count, spi, modifier);

  int stateSize = state.size();
  vector<Type> concatenated_state_vector((numIterations+1)*stateSize);
  vector<Type> concatenated_ratemat_nonzeros((numIterations+1)*updateidx.size());

  // Add initial state vector and non-zero element of the rate matrix into
  // corresponding vectors prefixed with "concatenated_".
  concatenated_state_vector.block(0, 0, stateSize, 1) = state;

  for (int j=0; j<updateidx.size(); j++) {
    int idx = updateidx[j] - 1;
    int row = from[idx] - 1;
    int col = to[idx] - 1;
    concatenated_ratemat_nonzeros[j] = ratemat.coeff(row,col);
  }

  // Calculate jacobian
  update_state_functor<Type> f(ratemat, par_accum_indices, do_hazard);
  matrix<Type> j = autodiff::jacobian(f, state);
  REPORT(j);

  int nextBreak = 0;
  int start = 0;
  for (int i=0; i<numIterations; i++) {

    state = do_step(state, ratemat, par_accum_indices, do_hazard);
    sp.block(0, 0, stateSize, 1) = state;

    // // update sp (state+params) and rate matrix
    // if (nextBreak<breaks.size() && i==(breaks[nextBreak])) {
    //     for (int j=start; j<start+count_of_tv_at_breaks[nextBreak]; j++) {
    //         //if (tv_orig[j])
    //         if (tv_update_method[j] & 0b00)
    //           // new value = original times another parameter
    //           sp[tv_spi[j]-1] = sp_orig[tv_spi[j]-1]*sp[tv_mult_spi[j]-1];
    //         else if ((tv_update_method[j] & 0b10) | (tv_update_method[j] & 0b11))
    //           // new value = previous times another parameter
    //           sp[tv_spi[j]-1] *= sp[tv_mult_spi[j]-1];
    //         else if (tv_update_method[j] & 0b01)
    //           // new value = another parameter
    //           sp[tv_spi[j]-1] = sp[tv_mult_spi[j]-1];
    //     }
    //
    //     start += count_of_tv_at_breaks[nextBreak];
    //     nextBreak++;
    // }

    // update sp (state+params) and rate matrix
    if (nextBreak<breaks.size() && i==(breaks[nextBreak])) {
      for (int j=start; j<start+count_of_tv_at_breaks[nextBreak]; j++) {
        if (tv_orig[j])
          sp[tv_spi[j]-1] = sp_orig[tv_spi[j]-1]*tv_mult[j];
        else
          sp[tv_spi[j]-1] *= tv_mult[j];
      }

      start += count_of_tv_at_breaks[nextBreak];
      nextBreak++;
    }

    // loop over a set of state variables, and update them with
    // expressions of other state variables and parameters

    update_sum_in_sp(sp, sumidx, sumcount, summandidx);
    //for (int j=0; j<sumidx.size(); j++) {
    //  std::cout << j << " idx = " << sumidx[j] << " sum = " << sp[sumidx[j]-1] << std::endl;
    //}
    //std::cout << "sp_110 = " << sp[110] << std::endl;
    update_ratemat(&ratemat, sp, from, to, count_integral, spi, modifier, updateidx);

    // concatenate state vectors at each time step so they can be returned
    concatenated_state_vector.block((i+1)*stateSize, 0, stateSize, 1) = state;

    // concatenate changing rate matrix elements at each time step so they can be returned
    int offset = (i+1)*updateidx.size();
    for (int j=0; j<updateidx.size(); j++) {
      int idx = updateidx[j] - 1;
      int row = from[idx] - 1;
      int col = to[idx] - 1;
      concatenated_ratemat_nonzeros[offset+j] = ratemat.coeff(row,col);
    }
  }

  //std::cout << concatenated_ratemat_nonzeros << std::endl;
  //std::cout << ratemat << std::endl;

  REPORT(ratemat);
  REPORT(concatenated_state_vector);
  REPORT(concatenated_ratemat_nonzeros);

  return state.sum();
}

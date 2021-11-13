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
// Helper function OutFlow
template<class Type>
vector<Type> OutFlow(
    const Eigen::SparseMatrix<Type>& mat,
    const vector<int>& outflow_row_count,
    const vector<int>& outflow_col_count,
    const vector<int>& outflow_rows,
    const vector<int>& outflow_cols)
{
  vector<Type> result(mat.rows());

  // Initialize result with zeros
  for (int i=0; i<result.size(); i++)
    result(i) = 0.0;

  int startRow = 0;
  for (int k=0; k<outflow_row_count.size(); k++) { // groups
    for (int i=startRow; i<startRow+outflow_row_count(k); i++) { // rows in a group
       int row = outflow_rows[i] - 1;

       int startCol = 0;
       for (int j=startCol; j<startCol+outflow_col_count(k); j++) { // cols in a row
         int col = outflow_cols[j] - 1;
         result[row] += mat.coeff(row, col);
       }
       startCol += outflow_cols(k);
    }
    startRow += outflow_row_count(k);
  }

  //std::cout << "++++++++++++++++++++++++++++" << std::endl;
  //std::cout << mat << std::endl;
  //std::cout << "----------------------------" << std::endl;
  //std::cout << result << std::endl;

  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Helper function Norm
template<class Type>
Type Norm(
    const vector<Type>& vec
)
{
  Type w = 0.0;
  for (int j=0; j<vec.size(); j++)
     w += vec(j)*vec(j);

  return sqrt(w);
}

///////////////////////////////////////////////////////////////////////////////
// Helper function CalcEigenVector
template<class Type>
vector<Type> CalcEigenVector(
    const matrix<Type>& jacobian,
    const vector<Type>& state,
    int iterations = 800,
    Type tolerance = 0.001)
{
  int n = state.size();

  //std::cout<< "n = " << n << std::endl;
  //std::cout<< "jacob = " << jacobian.rows() << ", " << jacobian.cols() << std::endl;
  //std::cout<< "state = " << state.rows() << ", " << state.cols() << std::endl;

  // Remove first and last two rows and columns from jacobian matrix and state vector
  matrix<Type> mat = jacobian.block(1, 1, n-3, n-3);
  vector<Type> vec = state.block(1, 0, n-3, 1);
  vector<Type> prevec(1);

  //std::cout<< "mat = " << mat << std::endl;
  //std::cout<< "vec 1 = " << vec << std::endl;

  int i;
  vector<Type> diff;

  for (i=0; i<iterations; i++) {
    vec = mat*vec;
    vec /= Norm(vec);

    if (i%50==0) {
      //std::cout << "======================= " << i << std::endl;
      if (prevec.size() != vec.size()) {
        prevec = vec;
      }
      else {
        diff = vec-prevec;

        if (Norm(diff) < tolerance) {
          //std::cout<< "diff = " << diff << std::endl;
          break;
        }
        else
          prevec = vec;
      }
    }
  }
  //std::cout << "==== Stop iteration at " << i << std::endl;
  //        std::cout<< "====pre vec = " << prevec << std::endl;
  //        std::cout<< "====cur vec (principla eigenvector) = " << vec << std::endl;
  //        std::cout<< "====diff = " << diff << std::endl;

  return vec;
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
//    vector<int> par_accum_indices,
    vector<int> outflow_row_count,
    vector<int> outflow_col_count,
    vector<int> outflow_rows,
    vector<int> outflow_cols,
    int do_hazard
)
{
  // Calculate flow matrix
  Eigen::SparseMatrix<Type> flows = calc_flowmat(ratemat, state, do_hazard);
  vector<Type> inflow = colSums(flows);
  //remove_cols(flows, par_accum_indices);
  //vector<Type> outflow = rowSums(flows); // remove some columns before doing so
  vector<Type> outflow = OutFlow(flows, outflow_row_count, outflow_col_count, outflow_rows, outflow_cols);
  state = state - outflow + inflow;
  return state;
}

///////////////////////////////////////////////////////////////////////////////
// Define Functor for jacobian

template <class Type>
struct update_state_functor{
  // Data members
  vector<Type> params_;
  vector<int> from_;
  vector<int> to_;
  vector<int> count_;
  vector<int> spi_;
  vector<int> modifier_;
  vector<int> sumidx_;
  vector<int> sumcount_;
  vector<int> summandidx_;
  // vector<int> par_accum_indices_;
  vector<int> linearized_outflow_row_count_;
  vector<int> linearized_outflow_col_count_;
  vector<int> linearized_outflow_rows_;
  vector<int> linearized_outflow_cols_;
  int do_hazard_;

  // Constructor
  update_state_functor(
    vector<Type> params,
    vector<int> from,
    vector<int> to,
    vector<int> count,
    vector<int> spi,
    vector<int> modifier,
    vector<int> sumidx,
    vector<int> sumcount,
    vector<int> summandidx,
    // vector<int> par_accum_indices,
    vector<int> linearized_outflow_row_count,
    vector<int> linearized_outflow_col_count,
    vector<int> linearized_outflow_rows,
    vector<int> linearized_outflow_cols,
    int do_hazard) : params_(params), from_(from), to_(to), count_(count),
                     spi_(spi), modifier_(modifier), sumidx_(sumidx), sumcount_(sumcount),
                     summandidx_(summandidx), // par_accum_indices_(par_accum_indices),
                     linearized_outflow_row_count_(linearized_outflow_row_count),
                     linearized_outflow_col_count_(linearized_outflow_col_count),
                     linearized_outflow_rows_(linearized_outflow_rows),
                     linearized_outflow_cols_(linearized_outflow_cols),
                     do_hazard_(do_hazard)
  {
  }

  // The function itself
  template <typename T>
  vector<T> operator()(vector<T> state_)
  {
    // Convert params_ from Type to T
    int n = params_.size();
    vector<T> params(n);
    for (int i=0; i<n; i++)
       params[i] = (T) params_[i];

    // Concatenate state and params
    vector<T> sp(state_.size()+params.size()+sumidx_.size());
    vector<T> place_holder(sumidx_.size());
    sp << state_, params, place_holder;

    update_sum_in_sp(sp, sumidx_, sumcount_, summandidx_);

    // We've got everything we need, lets do the job ...
    Eigen::SparseMatrix<T> ratemat = make_ratemat(state_.size(), sp, from_, to_, count_, spi_, modifier_);



    // 1 convert from Type to T
    //Eigen::SparseMatrix<T> ratemat;
    //ratemat = ratemat_.template cast<T>();

    // 2 do all the calculations in T
    vector<T> updated_state = do_step(state_, ratemat, // par_accum_indices_,
                                      linearized_outflow_row_count_, linearized_outflow_col_count_,
                                      linearized_outflow_rows_, linearized_outflow_cols_,
                                      do_hazard_);
    //std::cout << "updated_state = " << updated_state << std::endl;

    return (updated_state);
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
  //DATA_IVECTOR(par_accum_indices);

  DATA_IVECTOR(linearized_outflow_row_count);
  DATA_IVECTOR(linearized_outflow_col_count);
  DATA_IVECTOR(linearized_outflow_rows);
  DATA_IVECTOR(linearized_outflow_cols);

  DATA_IVECTOR(outflow_row_count);
  DATA_IVECTOR(outflow_col_count);
  DATA_IVECTOR(outflow_rows);
  DATA_IVECTOR(outflow_cols);

  DATA_VECTOR(lin_param_vals);
  DATA_IVECTOR(lin_param_count);
  DATA_IVECTOR(lin_param_idx);

  DATA_IVECTOR(df_state_par_idx);
  DATA_IVECTOR(df_state_count);
  DATA_IVECTOR(df_state_idx);

  DATA_IVECTOR(im_all_drop_eigen_idx);
  DATA_IVECTOR(im_eigen_drop_infected_idx);
  DATA_IVECTOR(im_all_to_infected_idx);
  DATA_IVECTOR(im_susceptible_idx);

  DATA_IVECTOR(ip_total_idx);
  DATA_IVECTOR(ip_infected_idx);

  DATA_IVECTOR(sumidx);
  DATA_IVECTOR(sumcount);
  DATA_IVECTOR(summandidx);

  DATA_INTEGER(numIterations);
  DATA_INTEGER(do_hazard);

  PARAMETER_VECTOR(params);

  // std::cout << "here in the objective function...";
  // state = make_state(params);

  std::cout << "linearized_outflow_row_count = " << linearized_outflow_row_count << std::endl;
  std::cout << "linearized_outflow_col_count = " << linearized_outflow_col_count << std::endl;
  std::cout << "linearized_outflow_rows = " << linearized_outflow_rows << std::endl;
  std::cout << "linearized_outflow_cols = " << linearized_outflow_cols << std::endl;

  std::cout << "outflow_row_count = " << outflow_row_count << std::endl;
  std::cout << "outflow_col_count = " << outflow_col_count << std::endl;
  std::cout << "outflow_rows = " << outflow_rows << std::endl;
  std::cout << "outflow_cols = " << outflow_cols << std::endl;

  std::cout << "lin_param_vals = " << lin_param_vals << std::endl;
  std::cout << "lin_param_count = " << lin_param_count << std::endl;
  std::cout << "lin_param_idx = " << lin_param_idx << std::endl;

  std::cout << "df_state_par_idx = " << df_state_par_idx << std::endl;
  std::cout << "df_state_count = " << df_state_count << std::endl;
  std::cout << "df_state_idx = " << df_state_idx << std::endl;

  std::cout << "im_all_drop_eigen_idx = " << im_all_drop_eigen_idx << std::endl;
  std::cout << "im_eigen_drop_infected_idx = " << im_eigen_drop_infected_idx << std::endl;
  std::cout << "im_all_to_infected_idx = " << im_all_to_infected_idx << std::endl;
  std::cout << "im_susceptible_idx = " << im_susceptible_idx << std::endl;

  std::cout << "ip_total_idx = " << ip_total_idx << std::endl;
  std::cout << "ip_infected_idx = " << ip_infected_idx << std::endl;

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
  //std::cout << ratemat << std::endl;
  //std::cout << par_accum_indices << std::endl;
  //std::cout << do_hazard << std::endl;
  //std::cout << state << std::endl;

  //update_state_functor<Type> f(ratemat, par_accum_indices, do_hazard);
  update_state_functor<Type> f(params, from, to, count, spi, modifier,
                               sumidx, sumcount, summandidx, // par_accum_indices,
                               linearized_outflow_row_count, linearized_outflow_col_count,
                               linearized_outflow_rows, linearized_outflow_cols, do_hazard);

  matrix<Type> j = autodiff::jacobian(f, state);

  //j = matrix<Type>::Random(3,3);
  vector<Type> eigenvec = CalcEigenVector(j, state, 5000);

//  Eigen::EigenSolver<Eigen::MatrixXd> es;
//  Eigen::MatrixXd A = Eigen::MatrixXd::Random(4,4);
//  es.compute(A, false);
//  std::cout << "The eigenvalues of A are: " << es.eigenvalues().transpose() << std::endl;


  REPORT(j);
  REPORT(eigenvec);

  int nextBreak = 0;
  int start = 0;
  for (int i=0; i<numIterations; i++) {

    state = do_step(state, ratemat, // par_accum_indices,
                    outflow_row_count, outflow_col_count,
                    outflow_rows, outflow_cols,
                    do_hazard);
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

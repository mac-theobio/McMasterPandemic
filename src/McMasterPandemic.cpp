// This version implements spec 0.1.1 in https://canmod.net/misc/flex_specs

// search for printouts that are not commented out: ^\s*std::cout
// search for printouts that are commented out: ^\s*//\s*std::cout

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/time.h>
#include <TMB.hpp>

///////////////////////////////////////////////////////////////////////////////
// Status of TMB calculations.
// If it is 0, then succeed. Otherwise, it encodes various causes for failure
//          1: doesn't not converge in CalcEigenVector;
//          2: eigen vector is all zeros in CalcEigenVector;
//          3: mixed signs in eigen vector in CalcEigenVector;

//static int tmb_status = 0;

///////////////////////////////////////////////////////////////////////////////
// Helper function round
// template<class Type>
// vector<int> round_vec(const vector<Type> vec_in)
// {
//   vector<int> vec_rounded(vec_in.size());
//   for (int i=0; i < vec_in.size(); i++)
//     vec_rounded[i] = (int) (vec_in[i] + 0.5);
//
//   return vec_rounded;
// }

///////////////////////////////////////////////////////////////////////////////
// Helper function rowSums -- don't think this will work
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
  int startCol = 0;
  for (int k=0; k<outflow_row_count.size(); k++) { // groups
    for (int i=startRow; i<startRow+outflow_row_count(k); i++) { // rows in a group
       int row = outflow_rows[i] - 1;
       for (int j=startCol; j<startCol+outflow_col_count(k); j++) { // cols in a row
         int col = outflow_cols[j] - 1;
         result[row] += mat.coeff(row, col);
       }
    }
    startCol += outflow_col_count(k);
    startRow += outflow_row_count(k);
  }
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
    int iterations = 8000,
    Type tolerance = 0.000001)
{
  //if (iterations<100)
  //  iterations = 100;	// this is the minimum

  // Remove first and last two rows and columns from jacobian matrix and state vector
  matrix<Type> mat = jacobian;
  vector<Type> vec = state;
  vector<Type> prevec(1);

  int i;
  vector<Type> diff;

  for (i=0; i<iterations; i++) {
    vec = mat*vec;
    vec /= Norm(vec);

    if (i%50==0) {

      if (prevec.size() != vec.size()) {
        prevec = vec;
      }
      else {
        diff = vec-prevec;

        // FIXME: parameter dependent branching??
        //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
        if (Norm(diff) < tolerance) {

          break;
        }
        else
          prevec = vec;
      }
    }
  }

  if (i==iterations) {
    //tmb_status = 1; 	// doesn't not converge
    return vec;
  }

  // check if the signs are the same
  //for (i=0; i<vec.size(); i++)
  //  if (vec[i]!=0) break;

  //if (i==vec.size()) {
  //  tmb_status = 2; 	// eigen vector is all zeros
  //  return vec;
  //}

  //for (int j=i+1; j<vec.size(); j++)
  //  if (vec[j-1]*vec[j]<0) {
  //    tmb_status = 3;     // mixed signs in eigen vector
  //    return vec;
  //  }

  // flip the sign
  //if (vec[i]<00)
  //  vec = -vec;

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
  //const double haz_eps_doub = 1e-12;
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
        //x = (1 - exp(-x / haz_eps_doub)) / (x + haz_eps_doub);
        if (x > 1e-12) {
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
  //const double haz_eps_doub = 1e-12;
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
        // FIXME: parameter dependent branching??
        //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
        //x = (1 - exp(-x / haz_eps_doub)) / (x + haz_eps_doub);
        if (x > 1e-12) {
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
    int do_hazard,
    int do_approx_hazard)
{
  if (!do_hazard) {
    return col_multiply(mat, vec);
  } else {
    const double haz_eps_doub = 1e-12;
    vector<Type> r = rowSums(mat);
    vector<Type> rho = exp(-r);
    vector<Type> s_tilde(vec.size());
    if (do_approx_hazard) {
      for (int i=0; i<vec.size(); i++) {
        s_tilde[i] = vec[i] * ((1 - exp(-r[i] / haz_eps_doub)) / (r[i] + haz_eps_doub));
      }
    } else {
      for (int i=0; i<vec.size(); i++) {
        // FIXME: parameter dependent branching??
        //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
        if (r[i]==0) // shall it be something like "abs(r[i])<0.00001" ?
         s_tilde[i] = 0.0;
        else
         s_tilde[i] = vec[i]/r[i];
      }
    }

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
    int do_hazard,
    int do_approx_hazard
)
{
  Eigen::SparseMatrix<Type> flows = calc_flowmat(ratemat, state, do_hazard, do_approx_hazard);
  vector<Type> inflow = colSums(flows);
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
  int do_approx_hazard_;

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
    int do_hazard,
    int do_approx_hazard) : params_(params), from_(from), to_(to), count_(count),
                     spi_(spi), modifier_(modifier), sumidx_(sumidx), sumcount_(sumcount),
                     summandidx_(summandidx), // par_accum_indices_(par_accum_indices),
                     linearized_outflow_row_count_(linearized_outflow_row_count),
                     linearized_outflow_col_count_(linearized_outflow_col_count),
                     linearized_outflow_rows_(linearized_outflow_rows),
                     linearized_outflow_cols_(linearized_outflow_cols),
                     do_hazard_(do_hazard),
                     do_approx_hazard_(do_approx_hazard)
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

    // do all the calculations in T
    vector<T> updated_state = do_step(state_, ratemat,
                                      linearized_outflow_row_count_, linearized_outflow_col_count_,
                                      linearized_outflow_rows_, linearized_outflow_cols_,
                                      do_hazard_, do_approx_hazard_);

    return (updated_state);
  }
};


///////////////////////////////////////////////////////////////////////////////
// We follow the steps specified in https://canmod.net/misc/flex_specs#v0.1.1
template<class Type>
vector<Type> make_state(
  const vector<Type>& params,
  int n_states,
  const vector<int>& lin_param_count,
  const vector<int>& lin_param_idx,
  const vector<Type>& lin_param_vals,
  const vector<int>& df_state_count,
  const vector<int>& df_state_idx,
  const vector<int>& df_state_par_idx,
  const vector<int>& from,
  const vector<int>& to,
  const vector<int>& count,
  const vector<int>& spi,
  const vector<int>& modifier,
  const vector<int>& sumidx,
  const vector<int>& sumcount,
  const vector<int>& summandidx,
  const vector<int>& linearized_outflow_row_count,
  const vector<int>& linearized_outflow_col_count,
  const vector<int>& linearized_outflow_rows,
  const vector<int>& linearized_outflow_cols,
  int do_hazard,
  int do_approx_hazard,
  const vector<int>& im_all_drop_eigen_idx,
  const vector<int>& im_eigen_drop_infected_idx,
  const vector<int>& im_all_to_infected_idx,
  const vector<int>& im_susceptible_idx,
  int ip_total_idx,
  int ip_infected_idx,
  int max_iters_eig_pow_meth,
  Type tol_eig_pow_meth
)
{

  int n_make_state_steps = 100;

  // 1 -- Initialize two state vectors,
  //      one for full model and one for linearized model
  vector<Type> state(n_states);
  vector<Type> lin_state(n_states);
  state = 0;
  lin_state = 0;

  // 2 -- Copy the parameters so they can be modified
  //      for the linearized model
  vector<Type> lin_params(params);

  // 3 -- Replace some elements of parameters for the
  //      linearized model
  int start = 0;
  for (int i=0; i<lin_param_count.size(); i++) {
    for (int j=start; j<start+lin_param_count[i]; j++) {
      lin_params[lin_param_idx[j]-1] = lin_param_vals[i];
    }
    start += lin_param_count[i];
  }

  // 4 -- Replace some elements of state for the
  //      linearized model
  start = 0;
  for (int i=0; i<df_state_count.size(); i++) {
    for (int j=start; j<start+df_state_count[i]; j++) {
      lin_state[df_state_idx[j]-1] = lin_params[df_state_par_idx[i]-1];
    }
    start += df_state_count[i];
  }

  // 5 -- Compute the Jacobian for the linearized model
  update_state_functor<Type> f(lin_params, from, to, count, spi, modifier,
                               sumidx, sumcount, summandidx,
                               linearized_outflow_row_count, linearized_outflow_col_count,
                               linearized_outflow_rows, linearized_outflow_cols,
                               do_hazard, do_approx_hazard);

  matrix<Type> jacob = autodiff::jacobian(f, lin_state);

  int nRows = jacob.rows();
  int nCols = jacob.cols();

  // 6 -- Remove rows and columns from the Jacobian
  //      and the state for the linearized model
  // Make a copy and append one nRows at the end.
  int n = im_all_drop_eigen_idx.size();
  vector<int> tmp_im_all_drop_eigen_idx(n+1);
  for (int i=0; i<n; i++)
    tmp_im_all_drop_eigen_idx[i] = im_all_drop_eigen_idx[i] - 1;	// 1-based indexing to 0-based indexing
  tmp_im_all_drop_eigen_idx[n] = nRows;		// Adding this at the end makes condensing step below be able to complete in one for loop

  //tmp_im_all_drop_eigen_idx[0] = 5;
  //tmp_im_all_drop_eigen_idx[1] = 0;

  n++;

  matrix<Type> trimmed_jacob(nRows-n+1, nCols-n+1);
  vector<Type> trimmed_lin_state(lin_state.size()-n+1);

  if (n>1) {
    // Sort the drop-out indices
    for (int i=n-1; i>0; i--)
      for (int j=0; j<i; j++)
        if (tmp_im_all_drop_eigen_idx[j]>tmp_im_all_drop_eigen_idx[j+1]) { // swap
          int t = tmp_im_all_drop_eigen_idx[j];
          tmp_im_all_drop_eigen_idx[j] = tmp_im_all_drop_eigen_idx[j+1];
          tmp_im_all_drop_eigen_idx[j+1] = t;
        }

    // Condense rows/columns by moving unremoved ones to left/up
    int empty_idx = tmp_im_all_drop_eigen_idx[0];
    for (int i=1; i<n; i++)
      for (int j=tmp_im_all_drop_eigen_idx[i-1]+1; j<tmp_im_all_drop_eigen_idx[i]; j++) {
        jacob.block(empty_idx, 0, 1, nCols) = jacob.block(j, 0, 1, nCols);
        jacob.block(0, empty_idx, nCols, 1) = jacob.block(0, j, nCols, 1);

        lin_state[empty_idx] = lin_state[j];

        empty_idx++;
      }

    // Copy the upper-left sub-matrix
    trimmed_jacob = jacob.block(0, 0, nRows-n+1, nCols-n+1);
    trimmed_lin_state = lin_state.block(0, 0, trimmed_lin_state.size(), 1);
  }

  // 7 -- Compute eigenvector
  vector<Type> eigenvec = CalcEigenVector(trimmed_jacob, trimmed_lin_state, max_iters_eig_pow_meth, tol_eig_pow_meth);

  // FIXME: parameter dependent branching??
  //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
  //if (tmb_status) return state; // There is an error in the computation so far

  // 8 -- Remove elements of `eigvec` to create `eig_infected`
  n = im_eigen_drop_infected_idx.size();
  vector<int> tmp_im_eigen_drop_infected_idx(n+1);
  for (int i=0; i<n; i++)
    tmp_im_eigen_drop_infected_idx[i] = im_eigen_drop_infected_idx[i] - 1;        // 1-based indexing to 0-based indexing
  tmp_im_eigen_drop_infected_idx[n] = eigenvec.size();         // Adding this at the end makes condensing step below be able to complete in one for loop

  n++;

  vector<Type> eig_infected(eigenvec.size()-n+1);

  if (n>1) {
    // Sort the drop-out indices
    for (int i=n-1; i>0; i--)
      for (int j=0; j<i; j++)
        if (tmp_im_eigen_drop_infected_idx[j]>tmp_im_eigen_drop_infected_idx[j+1]) { // swap
          int t = tmp_im_eigen_drop_infected_idx[j];
          tmp_im_eigen_drop_infected_idx[j] = tmp_im_eigen_drop_infected_idx[j+1];
          tmp_im_eigen_drop_infected_idx[j+1] = t;
        }


    // Condense rows/columns by moving unremoved ones to left/up
    int empty_idx = tmp_im_eigen_drop_infected_idx[0];
    for (int i=1; i<n; i++)
      for (int j=tmp_im_eigen_drop_infected_idx[i-1]+1; j<tmp_im_eigen_drop_infected_idx[i]; j++) {
        eigenvec[empty_idx] = eigenvec[j];

        empty_idx++;
      }

    eig_infected = eigenvec.block(0, 0, eig_infected.size(), 1);
  }

  // 9 -- Normalize `eig_infected` to have elements that sum to one
  eig_infected /= eig_infected.sum();

  // 10 -- distribute infected individuals among compartments in the initial state vector
  for (int i=0; i<im_all_to_infected_idx.size(); i++)
    state[im_all_to_infected_idx[i]-1] = eig_infected[i] * params[ip_infected_idx-1];

  // 11 -- distribute susceptible individuals among compartments in the initial state vector
  for (int i=0; i<im_susceptible_idx.size(); i++)
    state[im_susceptible_idx[i] - 1] = (1.0/im_susceptible_idx.size()) * (params[ip_total_idx-1] - params[ip_infected_idx-1]);

  return state;
}

///////////////////////////////////////////////////////////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Joint negative log-likelihood (stub)
  //Type jnll= 0;

  // Get data and parameters from R
  DATA_VECTOR(state);
  DATA_IVECTOR(from);
  DATA_IVECTOR(to);
  DATA_IVECTOR(count);
  DATA_IVECTOR(spi);
  DATA_IVECTOR(modifier);
  DATA_IVECTOR(updateidx);
  DATA_IVECTOR(breaks);
  DATA_IVECTOR(count_of_tv_at_breaks);
  DATA_IVECTOR(tv_spi);
  DATA_IVECTOR(tv_spi_unique);
  DATA_IVECTOR(tv_orig);
  DATA_IVECTOR(tv_abs);

  DATA_IVECTOR(outflow_row_count);
  DATA_IVECTOR(outflow_col_count);
  DATA_IVECTOR(outflow_rows);
  DATA_IVECTOR(outflow_cols);

  DATA_INTEGER(do_make_state);
  DATA_INTEGER(max_iters_eig_pow_meth);
  DATA_SCALAR(tol_eig_pow_meth);
  // DATA_SCALAR(haz_eps);  // TODO: not being used, but should be

  DATA_IVECTOR(linearized_outflow_row_count);
  DATA_IVECTOR(linearized_outflow_col_count);
  DATA_IVECTOR(linearized_outflow_rows);
  DATA_IVECTOR(linearized_outflow_cols);

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

  DATA_INTEGER(ip_total_idx);
  DATA_INTEGER(ip_infected_idx);

  DATA_IVECTOR(sumidx);
  DATA_IVECTOR(sumcount);
  DATA_IVECTOR(summandidx);

  DATA_INTEGER(numIterations);
  DATA_INTEGER(do_hazard);
  DATA_INTEGER(do_approx_hazard);
  DATA_INTEGER(do_hazard_lin);
  DATA_INTEGER(do_approx_hazard_lin);

  // The order of these PARAMETER_VECTOR macros
  // is important because it defines the order with
  // which non-fixed parameters are passed to the
  // objective function on the R side. From ?MakeADFun:
  // "The order of the PARAMETER_ macros defines the order
  // of parameters in the final objective function"
  PARAMETER_VECTOR(params);
  PARAMETER_VECTOR(tv_mult);

  //REPORT(tmb_status);

  // make state vector from params vector
  if (do_make_state) {
    state = make_state(
      params,
      state.size(),
      lin_param_count,
      lin_param_idx,
      lin_param_vals,
      df_state_count,
      df_state_idx,
      df_state_par_idx,
      from,
      to,
      count,
      spi,
      modifier,
      sumidx,
      sumcount,
      summandidx,
      linearized_outflow_row_count,
      linearized_outflow_col_count,
      linearized_outflow_rows,
      linearized_outflow_cols,
      do_hazard_lin,
      do_approx_hazard_lin,
      im_all_drop_eigen_idx,
      im_eigen_drop_infected_idx,
      im_all_to_infected_idx,
      im_susceptible_idx,
      ip_total_idx,
      ip_infected_idx,
      max_iters_eig_pow_meth,
      tol_eig_pow_meth
    );

  }

  // FIXME: parameter dependent branching??
  //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
  //if (tmb_status) return 0; // There is an error in the computation so far

  // Initialize the vectors that contain state, parameters,
  // and sums of these quantities
  vector<Type> sp(state.size()+params.size()+sumidx.size());
  vector<Type> place_holder(sumidx.size());
  sp << state, params, place_holder;
  vector<Type> sp_orig(sp); // sp_orig does not contain sums
  update_sum_in_sp(sp, sumidx, sumcount, summandidx);

  // Calculate integral of count
  vector<int> count_integral(count.size()+1);
  count_integral[0] = 0;
  for (int i=0; i<count.size(); i++)
    count_integral[i+1] = count_integral[i] + count[i];

  // We've got everything we need, lets do the job ...
  Eigen::SparseMatrix<Type> ratemat = make_ratemat(state.size(), sp, from, to, count, spi, modifier);

  // FIXME: parameter dependent branching??
  //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
  //if (tmb_status) return 0; // There is an error in the computation so far

  int stateSize = state.size();
  vector<Type> concatenated_state_vector((numIterations+1)*stateSize);
  vector<Type> concatenated_ratemat_nonzeros((numIterations+1)*updateidx.size());
  vector<Type> concatenated_time_varying_parameters((numIterations+1)*tv_spi_unique.size());

  // Add initial state vector and non-zero element of the rate matrix into
  // corresponding vectors prefixed with "concatenated_".
  concatenated_state_vector.block(0, 0, stateSize, 1) = state;

  for (int j=0; j<updateidx.size(); j++) {
    int idx = updateidx[j] - 1;
    int row = from[idx] - 1;
    int col = to[idx] - 1;
    concatenated_ratemat_nonzeros[j] = ratemat.coeff(row,col);
  }

  //if (tmb_status) {
  //  return 0; // There is an error in the computation so far
  //}

  int nextBreak = 0;
  int start = 0;
  for (int i=0; i<numIterations; i++) {

    state = do_step(state, ratemat,
                    outflow_row_count, outflow_col_count,
                    outflow_rows, outflow_cols,
                    do_hazard, do_approx_hazard);

    // FIXME: parameter dependent branching?? commenting out for now
    //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
    //if (tmb_status) {
    //  return 0; // There is an error in the computation so far
    //}

    sp.block(0, 0, stateSize, 1) = state;

    // update sp (state+params) and rate matrix
    if (nextBreak<breaks.size() && i==(breaks[nextBreak])) {
      for (int j=start; j<start+count_of_tv_at_breaks[nextBreak]; j++) {
        if (tv_abs[j]) { // type == 'abs'
          sp[tv_spi[j]-1] = tv_mult[j];
        }
        else if (tv_orig[j]) { // type == 'rel_orig'
          sp[tv_spi[j]-1] = sp_orig[tv_spi[j]-1]*tv_mult[j];
        }
        else { // type == 'rel_prev'
          sp[tv_spi[j]-1] *= tv_mult[j];
        }
      }

      start += count_of_tv_at_breaks[nextBreak];
      nextBreak++;
    }

    // loop over a set of state variables, and update them with
    // expressions of other state variables and parameters

    update_sum_in_sp(sp, sumidx, sumcount, summandidx);

    update_ratemat(&ratemat, sp, from, to, count_integral, spi, modifier, updateidx);

    // FIXME: parameter dependent branching??
    //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
    //if (tmb_status) {
    //  return 0; // There is an error in the computation so far
    //}

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

  REPORT(ratemat);
  REPORT(concatenated_state_vector);
  REPORT(concatenated_ratemat_nonzeros);

  return state.sum();
}

// This version implements spec 0.2.0 in https://canmod.net/misc/flex_specs

// grep for printouts that are not commented out: ^\s*std::cout
// grep for printouts that are commented out: ^\s*//\s*std::cout

// Spec 0.2.0
// Here is the order of columns that should go in the matrix:
// 1 State variables
// 2 Changing rate matrix elements (there is one such element in this model)
// 3 Sums
// 4 Factrs
// 5 Element-wise expressions involving the first four items
// 6 Lag-n differences of the first five items
// 7 Convolutions of the first five items

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h> 	// isnan() is defined
#include <sys/time.h>
#include <TMB.hpp>
#include <cppad/local/cond_exp.hpp>

const double EPSILON = 1.0e-10; // less than this value is considered as zero

///////////////////////////////////////////////////////////////////////////////
// To make debugging with gdb work, I define the following functions:
void print(vector<double> x) {
  std::cout << x << std::endl;
}
void print(vector<int> x) {
  std::cout << x << std::endl;
}

void print(matrix<double> x) {
  std::cout << x << std::endl;
}
void print(matrix<int> x) {
  std::cout << x << std::endl;
}

void print(array<double> x) {
  x.print();
}
void print(array<int> x) {
  x.print();
}

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
// Helper function
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
  matrix<Type> mat = jacobian;
  vector<Type> vec = state;
  vector<Type> prevec(1);

  int i;
  vector<Type> diff;

  for (i=0; i<iterations; i++) {
    vec = mat*vec;
    Type norm = Norm(vec);
    if (norm<EPSILON)
      Rf_error("Divided by zero in CalcEigenVector");

    vec /= norm;

    if (i%50==0) {
      if (prevec.size() != vec.size()) {
        prevec = vec;
      }
      else {
        diff = vec-prevec;
        //std::cout << "iteration = " << i << std::endl;
        //std::cout << "norm = " << Norm(diff) << std::endl;
        //std::cout << "tolerance = " << tolerance << std::endl;

        // FIXME: parameter dependent branching??
        //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
        // Type one = 1;
        // Type zero = 0;
        // std::cout << "cond exp = " << CppAD::CondExpLt(Norm(diff), tolerance, one, zero) << std::endl;

        if (Norm(diff) < tolerance) {
          break;
        }
        else
          prevec = vec;
      }
    }
  }

  if (i==iterations) {
    Rf_warning("The eigenvector did not converge.");
    return vec;
  }
  else { // It converges. Lets do a last check on negativity
    for (i=0; i<vec.size(); i++)
      if (vec[i]<0.0) {
        Rf_warning("Negative eigen vector");
        break;
      }
    return vec;
  }
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
      if (modifier[j] & 0b001) {
        x = 1-x;
      }
      else if (modifier[j] & 0b010) {
        //x = (1 - exp(-x / haz_eps_doub)) / (x + haz_eps_doub);
        if (x > 1e-12) {
           x = 1/x;
        }
      }
      prod *= x;
    }
    result.coeffRef(row,col) += prod;
    if (result.coeff(row,col)<0) {
      //std::cout << "negative rate mat " << result.coeff(row,col) << std::endl;
      //std::cout << "row " << row << std::endl;
      //std::cout << "col " << col << std::endl;
      Rf_warning("Negative rate matrix element");
    }

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
      if (modifier[j] & 0b001) {
        x = 1-x;
      }
      else if (modifier[j] & 0b010) {
        // FIXME: parameter dependent branching??
        //        https://github.com/kaskr/adcomp/wiki/Things-you-should-NOT-do-in-TMB
        //x = (1 - exp(-x / haz_eps_doub)) / (x + haz_eps_doub);
        if (x > 1e-12) {
          x = 1/x;
        }
      }
      prod *= x;
    }
    ratemat->coeffRef(row,col) += prod;
    if (ratemat->coeff(row,col)<0.0) {
      //std::cout << "negative rate mat " << ratemat->coeff(row,col) << std::endl;
      //std::cout << "row " << row << std::endl;
      //std::cout << "col " << col << std::endl;
      Rf_warning("Negative rate matrix element.");
    }
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

///////////////////////////////////////////////////////////////////////////////
// Helper function update_factr_in_sp
template<class Type>
void update_factr_in_sp(
    vector<Type>& sp,
    const vector<int>& factr_spi,
    const vector<int>& factr_count,
    const vector<int>& factr_spi_compute,
    const vector<int>& factr_modifier)
{

  int start = 0;
//^\s*std::cout << "start = " << start << std::endl;
  int n = factr_count.size();
//^\s*std::cout << "n = " << n << std::endl;
  for (int i=0; i<n; i++) {
    Type result = 0.0;
//^\s*std::cout << "result = " << result << std::endl;
//^\s*std::cout << "i = " << i << std::endl;
    int idx = factr_spi[i] - 1;
//^\s*std::cout << "idx = " << idx << std::endl;
    //std::cout << "sp before = " << std::endl;
    //std::cout << sp << std::endl;
    // sp(idx) = 0.0;
    Type prod = 1.0;
    for (int j=start; j<start+factr_count[i]; j++) {
//^\s*std::cout << "j = " << j << std::endl;
      //std::cout << "idx_compute[" << j << "] = " << factr_spi_compute[j]-1 << std::endl;
      Type x = sp(factr_spi_compute[j]-1);
//^\s*std::cout << "x = " << x << std::endl;
      //std::cout << "initial x = " << x << std::endl;
      //std::cout << "modifier = " << factr_modifier[j] << std::endl;
      if (factr_modifier[j] & 0b100) {
        //result.coeffRef(row,col) += prod;
        result += prod;
        prod = 1;
      }
      if (factr_modifier[j] & 0b001) {
        x = 1-x;
      } else if (factr_modifier[j] & 0b010) {
        //x = (1 - exp(-x / haz_eps_doub)) / (x + haz_eps_doub);
        if (x > 1e-12) {
          x = 1/x;
        }
      }
//^\s*std::cout << "x = " << x << std::endl;
      //std::cout << "x before prod = " << x << std::endl;
      prod *= x;
//^\s*std::cout << "prod = " << prod << std::endl;
    }
    //result.coeffRef(row,col) += prod;
    result +=  prod;
//^\s*std::cout << "result = " << result << std::endl;
    start += factr_count[i];
//^\s*std::cout << "start = " << start << std::endl;
    //std::cout << "result =  " << result << std::endl;
    //std::cout << "final sp = " << sp[idx] << std::endl;
    //std::cout << "sp after  = " << std::endl;
    //std::cout << sp << std::endl;
    sp[idx] = result;
//^\s*std::cout << "sp[idx] = " << sp[idx] << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////
// Helper function update_pow_in_sp
template<class Type>
void update_pow_in_sp(
    vector<Type>& sp,
    const vector<int>& powidx,
    const vector<int>& powarg1idx,
    const vector<int>& powarg2idx,
    const vector<int>& powconstidx)
{
  for (int j=0; j<powidx.size(); j++) {
    Type powarg1 = sp(powarg1idx[j]-1);
    Type powarg2 = sp(powarg2idx[j]-1);
    Type powconst = sp(powconstidx[j]-1);
    sp[powidx[j]-1] = powconst * pow(powarg1, powarg2);
  }
}

template<class Type>
vector<Type> do_step(
    vector<Type> state,
    Eigen::SparseMatrix<Type> ratemat,
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

  for (int i=0; i<state.size(); i++)
    if (state[i]<0.0) {
      Rf_warning("Negative state variable");
      break;
    }

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
  vector<int> factr_spi_;
  vector<int> factr_count_;
  vector<int> factr_spi_compute_;
  vector<int> factr_modifier_;
  vector<int> powidx_;
  vector<int> powarg1idx_;
  vector<int> powarg2idx_;
  vector<int> powconstidx_;
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
    vector<int> factr_spi,
    vector<int> factr_count,
    vector<int> factr_spi_compute,
    vector<int> factr_modifier,
    vector<int> powidx,
    vector<int> powarg1idx,
    vector<int> powarg2idx,
    vector<int> powconstidx,
    // vector<int> par_accum_indices,
    vector<int> linearized_outflow_row_count,
    vector<int> linearized_outflow_col_count,
    vector<int> linearized_outflow_rows,
    vector<int> linearized_outflow_cols,
    int do_hazard,
    int do_approx_hazard) : params_(params), from_(from), to_(to), count_(count),
                     spi_(spi), modifier_(modifier), sumidx_(sumidx), sumcount_(sumcount),
                     summandidx_(summandidx), // par_accum_indices_(par_accum_indices),
                     factr_spi_(factr_spi),
                     factr_count_(factr_count),
                     factr_spi_compute_(factr_spi_compute),
                     factr_modifier_(factr_modifier),
                     powidx_(powidx),
                     powarg1idx_(powarg1idx),
                     powarg2idx_(powarg2idx),
                     powconstidx_(powconstidx),
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
    vector<T> sp(state_.size()+params.size()+sumidx_.size()+factr_spi_.size()+powidx_.size());
    vector<T> place_holder(sumidx_.size());
    vector<T> place_holder_factr(factr_spi_.size());
    vector<T> place_holder_pow(powidx_.size());

    sp << state_, params, place_holder, place_holder_factr, place_holder_pow;

    update_sum_in_sp(sp, sumidx_, sumcount_, summandidx_);
    update_factr_in_sp(
      sp,
      factr_spi_,
      factr_count_,
      factr_spi_compute_,
      factr_modifier_
    );
    update_pow_in_sp(
      sp, powidx_, powarg1idx_, powarg2idx_, powconstidx_
    );

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
  const vector<int>& factr_spi,
  const vector<int>& factr_count,
  const vector<int>& factr_spi_compute,
  const vector<int>& factr_modifier,
  const vector<int>& powidx,
  const vector<int>& powarg1idx,
  const vector<int>& powarg2idx,
  const vector<int>& powconstidx,
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

  // int n_make_state_steps = 100;

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
                               factr_spi, factr_count,
                               factr_spi_compute,
                               factr_modifier,
                               powidx, powarg1idx, powarg2idx, powconstidx,
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
    trimmed_lin_state = lin_state.segment(0, trimmed_lin_state.size());
  }

  // 7 -- Compute eigenvector
  vector<Type> eigenvec = CalcEigenVector(trimmed_jacob, trimmed_lin_state, max_iters_eig_pow_meth, tol_eig_pow_meth);

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

    eig_infected = eigenvec.segment(0, eig_infected.size());
  }

  // 9 -- Normalize `eig_infected` to have elements that sum to one
  eig_infected /= eig_infected.sum();

  // 10 -- distribute infected individuals among compartments in the initial state vector
  for (int i=0; i<im_all_to_infected_idx.size(); i++) {
    state[im_all_to_infected_idx[i]-1] = eig_infected[i] * params[ip_infected_idx-1];
    if (state[im_all_to_infected_idx[i]-1]<0.0)
      Rf_warning("Negative state variable");
  }
  // 11 -- distribute susceptible individuals among compartments in the initial state vector
  for (int i=0; i<im_susceptible_idx.size(); i++) {
    state[im_susceptible_idx[i] - 1] = \
      (1.0/im_susceptible_idx.size()) * (params[ip_total_idx-1] - params[ip_infected_idx-1]);
    if (state[im_susceptible_idx[i] - 1]<0.0)
      Rf_warning("Negative state variable");
  }
  return state;
}

///////////////////////////////////////////////////////////////////////////////
// General Objective Function in spec 0.2.0
template<class Type>
class LossFunc {
public:
  int id;			// loss id
  std::vector<int> spi; 	// index to sp

  // This member function calculates and returns the loss
  Type run(const Type& observed, const Type& simulated, const vector<Type>& sp,
           int obs_do_sim_constraint, Type obs_sim_lower_bound) {
    Type var;
    Type lll = 0.0;
    // Type clamping_tolerance = 1e-12;
    Type clamped_simulated;
    Type param1, param2;

    // clamping
    if (obs_do_sim_constraint) // && simulated<obs_sim_lower_bound)
      clamped_simulated = simulated + obs_sim_lower_bound * \
                          (1.0/(1.0-(simulated-obs_sim_lower_bound)/obs_sim_lower_bound + \
                                ((simulated-obs_sim_lower_bound)*(simulated-obs_sim_lower_bound))/(obs_sim_lower_bound*obs_sim_lower_bound) \
                               ) \
                          );
    else
      clamped_simulated = simulated;

    switch (id) {
      case 0: // Negative Binomial Negative Log Likelihood
        // var = mu + mu^2/k
        // p.165: https://ms.mcmaster.ca/~bolker/emdbook/book.pdf
        //   var ~ variance
        //   mu ~ mean
        //   k ~ overdispersion parameter = sp[this->spi[0]]
        var = clamped_simulated + ((clamped_simulated*clamped_simulated) / sp[this->spi[0]]);

        // k = (mu^2)/(var - mu)

        //std::cout << "var = " << var << std::endl;
        if (clamped_simulated<=0.0)
          Rf_warning("negative simulation value being compared with data");

        if (var<=0.0)
          Rf_warning("negative binomial variance is negative");

        if (sp[this->spi[0]] > 1e10)
          Rf_warning("negative binomial dispersion is very large, >1e10");

        // when k is large this should tend to log(factorial(observed)),
        // but we have observed NaNs in this case. not sure how to smoothly
        // (i.e. differentiably) address this
        lll = -1.0 * dnbinom2(observed, clamped_simulated, var, 1);
        //std::cout << "obs = " << observed << " sim = " << clamped_simulated << " var = " << var << " disp = " << sp[this->spi[0]] << " loss = " << lll << std::endl;
        break;

      case 1: // Poisson
        lll = -1.0 * dpois(observed, clamped_simulated, 1);
        break;

      case 2: // Log Normal
        param1 = sp[this->spi[0]]; // sd
        lll = -1.0 * dnorm(log(observed), log(clamped_simulated), param1, 1);
        break;

      case 3: // Normal
        param1 = sp[this->spi[0]]; // sd
        lll = -1.0 * dnorm(observed, clamped_simulated, param1, 1);
        break;

      case 4: // Beta
        Rf_error("Beta distributed error is under construction");
        param1 = sp[this->spi[0]]; // shape1
        param2 = sp[this->spi[1]]; // shape2
        //lll = pbeta(observed, param1, param2);
        break;

      default:
        // tmb error
        Rf_error("Unrecognized loss id");
        return 0.0;
    }

    // If there is no error
    return lll;
  };

  Type sim(const Type& simulated, const vector<Type>& sp,
           int obs_do_sim_constraint, Type obs_sim_lower_bound) {
    Type var;
    Type sims_with_error = simulated;
    Type clamped_simulated;
    Type var_tolerance = 1e-8;
    Type param1;

    // clamping
    if (obs_do_sim_constraint) // && simulated<obs_sim_lower_bound)
      clamped_simulated = simulated + obs_sim_lower_bound * \
                          (1.0/(1.0-(simulated-obs_sim_lower_bound)/obs_sim_lower_bound + \
                                ((simulated-obs_sim_lower_bound)*(simulated-obs_sim_lower_bound))/(obs_sim_lower_bound*obs_sim_lower_bound) \
                               ) \
                          );
    else
      clamped_simulated = simulated;

    switch (id) {
      case 0: // Negative Binomial Negative Log Likelihood
        // var = mu + mu^2/k
        // p.165: https://ms.mcmaster.ca/~bolker/emdbook/book.pdf
        //   var ~ variance
        //   mu ~ mean
        //   k ~ overdispersion parameter = sp[this->spi[0]]

        var = clamped_simulated + ((clamped_simulated*clamped_simulated) / sp[this->spi[0]]);

        if (var < var_tolerance)
          // more numerically stable to just set the simulations
          // to the mean when the var is low
          sims_with_error = clamped_simulated;
        else
          sims_with_error = rnbinom2(clamped_simulated, var);

        break;

      case 1: // Poisson
        sims_with_error = rpois(clamped_simulated);
        break;

      case 2: // Log Normal
        param1 = sp[this->spi[0]]; // sd
        sims_with_error = exp(rnorm(log(clamped_simulated), param1));
        break;

      case 3: // Normal
        param1 = sp[this->spi[0]]; // sd
        sims_with_error = rnorm(clamped_simulated, param1);
        break;

      case 4: // Beta
        //sims_with_error = rpois(clamped_simulated);
        break;

      default:
        // if you can't find a valid error distribution,
        // just return the mean
        break;
    }

    return sims_with_error;
  }
};

///////////////////////////////////////////////////////////////////////////////
// Regularization
template<class Type>
Type Regularization(
  const vector<Type>& params, // either params or tv_mult
  const vector<int>& param_id,
  const vector<int>& reg_family_id,
  const vector<int>& count_reg_params,
  const vector<Type>& reg_params
)
{
  // Shall we adjust params a bit to make sure positive value is passed to the log function?
  // const Type EPSILON = 1.0e-8;

  int start = 0;
  Type result = 0.0;

  for (int j=0; j<param_id.size(); j++) {
    // Based on the regularization type
    switch (reg_family_id[j]) {
      case 1: // Flat
        break;

      case 2: // Normal
        result += -(dnorm(params[param_id[j]-1], reg_params[start], reg_params[start+1], 1));
        break;

      // case 3: // others ...

      default: // error
        Rf_error("Unrecognized regularization type");
        break;
    }

    start += count_reg_params[j];
  }

  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Inverse transform
template<class Type>
void InverseTransformParams(
  vector<Type>& params, // either params or tv_mult
  const vector<int>& param_id,
  const vector<int>& trans_id
)
{
  //std::cout << "InverseTransformParams ..." << std::endl;
  //std::cout << param_id << std::endl;
  //std::cout << trans_id << std::endl;

  for (int j=0; j<param_id.size(); j++) {
    // Transformation no matter it is flat or not

    switch (trans_id[j]) {
      case 1: // identity
        //transformed = params[param_id[j]-1];
        break;
      case 2: // inverse log
        params[param_id[j]-1] = exp(params[param_id[j]-1]);
        //std::cout << params[param_id[j]-1] << std::endl;
        break;
      case 3: // inverse log_10
        params[param_id[j]-1] = pow(10.0, params[param_id[j]-1]);
        break;
      case 4: // inverse logit
        params[param_id[j]-1] = invlogit(params[param_id[j]-1]); // 1.0/(1.0+exp(-params[param_id[j]-1]));
        break;
      case 5: // inverse cloglog
        params[param_id[j]-1] = 1.0 - exp(-exp(params[param_id[j]-1]));
        break;
      case 6: // inverse of inverse
        if (params[param_id[j]-1]<1.0e-8 || params[param_id[j]-1]>-1.0e-8) {
          Rf_error("zero in inverse trans");
        }
        params[param_id[j]-1] = 1.0/params[param_id[j]-1];
        break;
      default: // unrecognized trans_id. Treat it as identity
        Rf_error("Unrecognized parameter transformation type");
        break;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{

  //std::cout << " =========== inside TMB ===========" << std::endl;

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
  DATA_IVECTOR(tv_type_id);

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

  DATA_IVECTOR(factr_spi);
  DATA_IVECTOR(factr_count);
  DATA_IVECTOR(factr_spi_compute);
  DATA_IVECTOR(factr_modifier);

  DATA_IVECTOR(powidx);
  DATA_IVECTOR(powarg1idx);
  DATA_IVECTOR(powarg2idx);
  DATA_IVECTOR(powconstidx);

  DATA_INTEGER(numIterations);
  DATA_INTEGER(do_hazard);
  DATA_INTEGER(do_approx_hazard);
  DATA_INTEGER(do_hazard_lin);
  DATA_INTEGER(do_approx_hazard_lin);

  DATA_IVECTOR(sr_count); 	// condense_count
  DATA_IVECTOR(sri)		// condense_sri
  DATA_IVECTOR(sr_modifier); 	// condense_modifier

  DATA_IVECTOR(lag_diff_sri);
  DATA_IVECTOR(lag_diff_delay_n);

  DATA_IVECTOR(conv_sri);
  DATA_IVECTOR(conv_c_prop_idx);
  DATA_IVECTOR(conv_c_delay_cv_idx);
  DATA_IVECTOR(conv_c_delay_mean_idx);
  DATA_IVECTOR(conv_qmax);

  DATA_IVECTOR(obs_var_id);
  DATA_IVECTOR(obs_loss_id);
  DATA_IVECTOR(obs_loss_param_count);
  DATA_IVECTOR(obs_spi_loss_param);
  DATA_IVECTOR(obs_time_step);
  DATA_IVECTOR(obs_history_col_id);
  DATA_VECTOR(obs_value);
  DATA_INTEGER(obs_do_sim_constraint);
  DATA_SCALAR(obs_sim_lower_bound);

  DATA_IVECTOR(opt_param_id);
  DATA_IVECTOR(opt_trans_id);
  DATA_IVECTOR(opt_count_reg_params);
  DATA_VECTOR(opt_reg_params);
  DATA_IVECTOR(opt_reg_family_id);

  //std::cout << "opt_param_id = " << opt_param_id << std::endl;
  //std::cout << "opt_trans_id = " << opt_trans_id << std::endl;
  //std::cout << "opt_count_reg_params = " << opt_count_reg_params << std::endl;
  //std::cout << "opt_reg_params = " << opt_reg_params << std::endl;
  //std::cout << "opt_reg_family_id = " << opt_reg_family_id << std::endl;

  DATA_IVECTOR(opt_tv_param_id);
  DATA_IVECTOR(opt_tv_trans_id);
  DATA_IVECTOR(opt_tv_count_reg_params);
  DATA_VECTOR(opt_tv_reg_params);
  DATA_IVECTOR(opt_tv_reg_family_id);

  //std::cout << "opt_tv_param_id = " << opt_tv_param_id << std::endl;
  //std::cout << "opt_tv_trans_id = " << opt_tv_trans_id << std::endl;

  // used for testing convolution code only
  //vector<int> conv_qmax(1); // you need to comment out DATA_IVECTOR(conv_qmax);
  //conv_qmax(0) = 6;
  //numIterations = 35;

  // The order of these PARAMETER_VECTOR macros
  // is important because it defines the order with
  // which non-fixed parameters are passed to the
  // objective function on the R side. From ?MakeADFun:
  // "The order of the PARAMETER_ macros defines the order
  // of parameters in the final objective function"
  PARAMETER_VECTOR(params);
  PARAMETER_VECTOR(tv_mult);

  // Regularization
  Type regularization = Regularization(
                          params,
                          opt_param_id,
                          opt_reg_family_id,
                          opt_count_reg_params,
                          opt_reg_params
                        );

  regularization += Regularization(
                      tv_mult,
                      opt_tv_param_id,
                      opt_tv_reg_family_id,
                      opt_tv_count_reg_params,
                      opt_tv_reg_params
                    );

  InverseTransformParams(
    params,
    opt_param_id,
    opt_trans_id
  );
  //std::cout << "params: " << params << std::endl;

  InverseTransformParams(
    tv_mult,
    opt_tv_param_id,
    opt_tv_trans_id
  );
  //std::cout << "tv_mult: " << tv_mult << std::endl;

  // dump inverse-transformed parameters and time-varying
  // multipliers to the R side, as an easy way to
  // back-transform the parameters off of the scale used
  // for optimization
  REPORT(params);
  REPORT(tv_mult);

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
      factr_spi,
      factr_count,
      factr_spi_compute,
      factr_modifier,
      powidx,
      powarg1idx,
      powarg2idx,
      powconstidx,
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

  // numbers of different types of columns in the simulation_history
  //  -- spec v0.2.0
  int stateSize = state.size(); // state variables
  int tvElementsNum = updateidx.size(); // time-varying rate matrix elements
  int sumSize = sumidx.size(); // sums of parameters and state variables
  int factrSize = factr_spi.size(); // factrs
  int powSize = powidx.size(); // powers
  int extraExprNum = sr_count.size(); // additional expressions
  int lagNum = lag_diff_sri.size(); // lagged differences
  int convNum = conv_sri.size(); // convolutions

  // std::cout << "powSize: " << powSize << std::endl;

  matrix<Type> simulation_history(numIterations+1, \
    stateSize+tvElementsNum+sumSize+factrSize+powSize+extraExprNum+lagNum+convNum);
  simulation_history.setZero();

  simulation_history.block(0, 0, 1, stateSize) = state.transpose();

  // Initialize the vectors that contain state, parameters,
  // and sums of these quantities
  vector<Type> sp(state.size()+params.size()+sumSize+factrSize+powSize);
  vector<Type> place_holder(sumSize);
  vector<Type> place_holder_factr(factrSize);
  vector<Type> place_holder_pow(powSize);
  sp << state, params, place_holder, place_holder_factr, place_holder_pow;
  vector<Type> sp_orig(sp); // sp_orig does not contain sums

  update_sum_in_sp(sp, sumidx, sumcount, summandidx);
  if (sumSize>0)
    simulation_history.block(0, stateSize+tvElementsNum, 1, sumSize) = \
      sp.segment(stateSize+params.size(), sumSize).transpose();

  update_factr_in_sp(
    sp,
    factr_spi, factr_count,
    factr_spi_compute, factr_modifier);
  if (factrSize>0)
    simulation_history.block(0, stateSize+tvElementsNum+sumSize, 1, factrSize) = \
      sp.segment(stateSize+params.size()+sumSize, factrSize).transpose();

  update_pow_in_sp(sp, powidx, powarg1idx, powarg2idx, powconstidx);
  if (powSize>0)
    simulation_history.block(0, stateSize+tvElementsNum+sumSize+factrSize, 1, powSize) = \
      sp.segment(stateSize+params.size()+sumSize+factrSize, powSize).transpose();

  // Calculate integral of count
  vector<int> count_integral(count.size()+1);
  count_integral[0] = 0;
  for (int i=0; i<count.size(); i++)
    count_integral[i+1] = count_integral[i] + count[i];

  // We've got everything we need, lets do the job ...
  Eigen::SparseMatrix<Type> ratemat = make_ratemat(state.size(), sp, from, to, count, spi, modifier);

  vector<Type> concatenated_state_vector((numIterations+1)*stateSize);
  vector<Type> concatenated_ratemat_nonzeros((numIterations+1)*updateidx.size());

  // Add initial state vector and non-zero element of the rate matrix into
  // corresponding vectors prefixed with "concatenated_".
  concatenated_state_vector.segment(0, stateSize) = state;

  for (int j=0; j<updateidx.size(); j++) {
    int idx = updateidx[j] - 1;
    int row = from[idx] - 1;
    int col = to[idx] - 1;
    concatenated_ratemat_nonzeros[j] = ratemat.coeff(row,col);
    simulation_history(0, stateSize+j) = ratemat.coeff(row,col);
  }

  // Item #5 Element-wise sum of any variable of type 4 (similar to factr calculation)
  if (extraExprNum>0) {
    vector<int> sr_output_idx(extraExprNum);
    for (int k=0; k<extraExprNum; k++) {
      sr_output_idx(k) = stateSize+tvElementsNum+sumSize+factrSize+powSize+k+1; // 1-based indexing
    }

    vector<Type> sh_row = simulation_history.row(0).transpose();

    update_factr_in_sp(sh_row, //simulation_history.row(i+1),
                       sr_output_idx,
                       sr_count,
                       sri,
                       sr_modifier);

    simulation_history.row(0) = sh_row.transpose();
  }

  // Item #7 preparation --- calculating kappa so that we don't need to
  // repeatedly calculate it in the simulation
  vector<vector<Type> > kappa(conv_sri.size());
  for (int k=0; k<conv_sri.size(); k++) {
    if (conv_qmax[k]<2) continue; // 2 is the mininum

    //std::cout << "kappa initial len=" << kappa[k].size() << std::endl;

    Type c_prop = params(conv_c_prop_idx[k]-1);
    Type c_delay_cv   = params(conv_c_delay_cv_idx[k]-1);
    Type c_delay_mean = params(conv_c_delay_mean_idx[k]-1);

    //std::cout << "conv_c_delay_cv_idx[k]=" << conv_c_delay_cv_idx[k] << std::endl;
    //std::cout << "c_delay_cv=" << c_delay_cv << std::endl;

    Type shape = 1.0/(c_delay_cv*c_delay_cv);
    Type scale = c_delay_mean/shape;

    vector<Type> delta(conv_qmax[k]-1);

    //std::cout << "shape=" << shape << std::endl;
    //std::cout << "scale=" << scale << std::endl;

    Type pre_gamma = pgamma ((Type) 1.0, shape, scale);
    for (int q=1; q<conv_qmax[k]; q++) {
      //std::cout << pre_gamma << std::endl;
      Type cur_gamma = pgamma ((Type) (q+1), shape, scale);
      delta(q-1) = cur_gamma - pre_gamma;
      pre_gamma = cur_gamma;
    }
    //std::cout << "delta = " << delta << std::endl;
    //std::cout << "delta sum = " << delta.sum() << std::endl;
    //std::cout << "c_prop = " << c_prop << std::endl;
    kappa[k] = c_prop*delta/delta.sum();
    //std::cout << "kappa = " << kappa[k] << std::endl;
  }

  // Preparation for the General Objective Function in spec 0.2.0
  // Build a map of var_id to loss func
  int start = 0;
  std::map<int, LossFunc<Type> > varid2lossfunc;
  for (int i=0; i<obs_var_id.size(); i++) {
    LossFunc<Type> lf;
    lf.id = obs_loss_id[i] - 1;

    for (int j=0; j<obs_loss_param_count[i]; j++)
      lf.spi.push_back(obs_spi_loss_param[start+j]-1);
    start += obs_loss_param_count[i];

    varid2lossfunc[obs_var_id[i]-1] = lf;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Simulation loop through "numIterations" time steps
  int nextBreak = 0;
  start = 0;
  Type sum_of_loss = 0.0;
  int obs_start = 0;

  for (int i=0; i<numIterations; i++) {

    //std::cout << "sp:" << std::endl;
    //std::cout << sp << std::endl;

    //std::cout << "ratemat:" << std::endl;
    //std::cout << ratemat << std::endl;

    state = do_step(state, ratemat,
                    outflow_row_count, outflow_col_count,
                    outflow_rows, outflow_cols,
                    do_hazard, do_approx_hazard);

    sp.segment(0, stateSize) = state;
    simulation_history.block(i+1, 0, 1, stateSize) = state.transpose();

    // update sp (params)
    if (nextBreak<breaks.size() && i==(breaks[nextBreak])) {
      for (int j=start; j<start+count_of_tv_at_breaks[nextBreak]; j++) {
        // valid_tv_types = c(
        //   'abs',
        //   'rel_orig', 'rel_prev',
        //   'rel_orig_logit', 'rel_prev_logit'
        // )
        switch (tv_type_id[j]-1) {
          case 0: // abs
            sp[tv_spi[j]-1] = tv_mult[j];
            break;
          case 1: // rel_orig
            sp[tv_spi[j]-1] = sp_orig[tv_spi[j]-1]*tv_mult[j];
            break;
          case 2: // rel_prev
            sp[tv_spi[j]-1] *= tv_mult[j];
            break;
          case 3: // rel_orig_logit
            sp[tv_spi[j]-1] = invlogit(
              logit(sp_orig[tv_spi[j]-1]) + tv_mult[j]
            );
            break;
          case 4: // rel_prev_logit
            sp[tv_spi[j]-1] = invlogit(
              logit(sp[tv_spi[j]-1]) + tv_mult[j]
            );
            break;
          default:
            Rf_error("invalid time-variation specification");
        }
        //std::cout << "tv_spi = " << tv_spi[j] << std::endl;
        //std::cout << "tv_mult = " << tv_mult[j] << std::endl;
        // if (tv_abs[j]) { // type == 'abs'
        //   sp[tv_spi[j]-1] = tv_mult[j];
        // }
        // else if (tv_orig[j]) { // type == 'rel_orig'
        //   sp[tv_spi[j]-1] = sp_orig[tv_spi[j]-1]*tv_mult[j];
        // }
        // else { // type == 'rel_prev'
        //   sp[tv_spi[j]-1] *= tv_mult[j];
        // }
      }

      start += count_of_tv_at_breaks[nextBreak];
      nextBreak++;
    }

    // loop over a set of state variables, and update them with
    // expressions of other state variables and parameters

    update_sum_in_sp(sp, sumidx, sumcount, summandidx);

    if (sumSize>0)
      simulation_history.block(i+1, stateSize+tvElementsNum, 1, sumSize) = \
        sp.segment(stateSize+params.size(), sumSize).transpose();

    update_factr_in_sp(
      sp, factr_spi,
      factr_count,
      factr_spi_compute,
      factr_modifier);

    if (factrSize>0)
      simulation_history.block(i+1, stateSize+tvElementsNum+sumSize, 1, factrSize) = \
        sp.segment(stateSize+params.size()+sumSize, factrSize).transpose();

    update_pow_in_sp(sp, powidx, powarg1idx, powarg2idx, powconstidx);
    if (powSize>0)
      simulation_history.block(i+1, stateSize+tvElementsNum+sumSize+factrSize, 1, powSize) = \
        sp.segment(stateSize+params.size()+sumSize+factrSize, powSize).transpose();

    update_ratemat(
      &ratemat, sp, from, to, count_integral,
      spi, modifier, updateidx);

    // concatenate state vectors at each time step so they can be returned
    concatenated_state_vector.segment((i+1)*stateSize, stateSize) = state;

    // concatenate changing rate matrix elements at each time step so they can be returned
    int offset = (i+1)*updateidx.size();
    for (int j=0; j<updateidx.size(); j++) {
      int idx = updateidx[j] - 1;
      int row = from[idx] - 1;
      int col = to[idx] - 1;
      concatenated_ratemat_nonzeros[offset+j] = ratemat.coeff(row,col);
      simulation_history(i+1, stateSize+j) = ratemat.coeff(row,col);
    }

    // Item #5 Element-wise sum of any variable of type 4 (similar to factr calculation)
    if (extraExprNum>0) {
      vector<int> sr_output_idx(extraExprNum);
      for (int k=0; k<extraExprNum; k++) {
        sr_output_idx(k) = stateSize+tvElementsNum+sumSize+factrSize+powSize+k+1; // 1-based indexing
      }

      vector<Type> sh_row = simulation_history.row(i+1).transpose();

      update_factr_in_sp(sh_row, //simulation_history.row(i+1),
                         sr_output_idx,
                         sr_count,
                         sri,
                         sr_modifier);

      simulation_history.row(i+1) = sh_row.transpose();
    }

    // The above version calls update_factr_in_sp function to do item #5. It's good
    // to reuse the code. The downside is that we need to transpose vector twice
    // and two extra vector assignments.
    //
    // Below I copy and edit the body of update_factr_in_sp, which is more efficient
    // but with a lot of duplicated code.
    // Note: We need to sp->simulation_history's row.  The version below is not fully tested.
    /*
    int start = 0;
    int n = sr_count.size();
    for (int k=0; k<n; k++) {
      Type result = 0.0;
      Type prod = 1.0;
      for (int j=start; j<start+sr_count[k]; j++) {
        Type x = sp(sri[j]-1);
        if (sr_modifier[j] & 0b100) {
          result += prod;
          prod = 1;
        }
        if (sr_modifier[j] & 0b001) {
          x = 1-x;
        } else if (sr_modifier[j] & 0b010) {
          if (x > 1e-12) {
            x = 1/x;
          }
        }
        prod *= x;
      }
      result +=  prod;
      start += sr_count[k];
      simulation_history(i+1, stateSize+tvElementsNum+sumSize+factrSize+k) = result;
    }
    */

    // Item #6 Lag-n differences of any variables of type 1-5
    for (int k=0; k<lag_diff_sri.size(); k++) {
      if (i+1>=lag_diff_delay_n[k]) {
        int col = lag_diff_sri[k]-1;
        simulation_history(i+1, stateSize+tvElementsNum+sumSize+factrSize+powSize+extraExprNum+k) = \
          simulation_history(i+1, col) - simulation_history(i+1-lag_diff_delay_n[k], col);
      }
    }

    // Item #7 Convolutions of any variables of type 1-5 with a gamma-density kernel
    int index_to_item7 = stateSize + tvElementsNum + sumSize + factrSize + powSize + \
                         extraExprNum + lag_diff_sri.size();
    for (int k=0; k<conv_sri.size(); k++) {
      vector<Type> kernel = kappa[k];
      //std::cout << "THIS IS REAL KERNEL = " << kernel << std::endl;
      Type conv = 0.0;
      if (i>conv_qmax[k]-4) { // i+2>=qmax-1
        //std::cout << "========= i = " << i << std::endl;
        for (int j=0; j<conv_qmax[k]-1; j++) {
          //std::cout << "x = " << simulation_history(i+1-j, conv_sri[k]) << std::endl;
          //std::cout << "k = " << kernel(j) << std::endl;
          conv += simulation_history(i+1-j, conv_sri[k]-1) * kernel(j);
          //std::cout << "z = " << conv << std::endl;
        }
        simulation_history(i+1, index_to_item7+k) = conv;
      }
    }

    // Calculate the General Objective Function
    for (int k=obs_start; k<obs_time_step.size(); k++) {
      if (i+1==obs_time_step[k]-1) {
        //std::cout << "k= " << k << " [" << obs_time_step[k] << ", " << obs_history_col_id[k] << "]" << std::endl;
        Type x = obs_value[k];
        Type mu = simulation_history(obs_time_step[k]-1, obs_history_col_id[k]-1);
        sum_of_loss += varid2lossfunc[obs_history_col_id[k]-1].run(x, mu, sp, obs_do_sim_constraint, obs_sim_lower_bound);
        //std::cout << "loglik = " << sum_of_loss << std::endl;
      }
      else {
        obs_start = k;
        break;
      }
    }


  }

  for (int i=0; i<numIterations+1; i++) {
    SIMULATE {
      for (int k=0; k<obs_var_id.size(); k++) {
        simulation_history(i, obs_var_id[k] - 1) = varid2lossfunc[obs_var_id[k] - 1].sim(
          simulation_history(i, obs_var_id[k] - 1),
          sp,
          obs_do_sim_constraint,
          obs_sim_lower_bound
        );
      }
    }
  }

  //std::cout << "Loss = " << sum_of_loss << std::endl;

  // Regularization
  sum_of_loss += regularization;
  //std::cout << "regularization = " << sum_of_loss << std::endl;

  //std::cout << "simulation_history size= " << simulation_history.size() << std::endl;
  //std::cout << simulation_history << std::endl;


  REPORT(ratemat);
  REPORT(concatenated_state_vector);
  REPORT(concatenated_ratemat_nonzeros);
  REPORT(simulation_history);

  return sum_of_loss;
}
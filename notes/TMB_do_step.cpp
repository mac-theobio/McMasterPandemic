#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/time.h>
#include <TMB.hpp>

// You must "sudo apt install r-mathlib" 
//#define MATHLIB_STANDALONE 1
//#include "Rmath.h"

///////////////////////////////////////////////////////////////////////////////
// In refactoring McMasterPandemic R code with TMB, I will keep code structure
// unchanged and use the same names of variables/functions so that
// (1) It is easier to compare the R code with the C++ code and find difference
// (2) It is easier to update the C++ code whenever the R code get updated.
//////////////////////////////////////////////////////////////////////////////

//#define TMB_VERBOSE

// hard-code the names of state and parameters in fixed order
static std::vector<std::string> stateNames {"S", "E", "Ia", "Ip", "Im", "Is", "H", "H2", \
                                "ICUs", "ICUd", "D", "R", "X", "V"};
// Perhaps we don't need parameter names, but I just keep it for possible future use
static std::vector<std::string> paramNames {"beta0", "Ca", "Cp", "Cm", "Cs" \
                                "alpha", "sigma", "gamma_a", "gamma_m", "gamma_s", \
                                "gamma_p", "rho", "delta", "mu", "N" \
                                "E0", "nonhosp_mort", "iso_m", "iso_s", "phi1", \
                                "phi2", "psi1", "psi2", "psi3", "c_prop", \
                                "c_delay_mean", "c_delay_cv", "proc_disp", "zeta"};



///////////////////////////////////////////////////////////////////////////////
// Helper function exclude_states
// Remove every element of vector "b" from vector "a"
std::vector<std::string> exclude_states(
    const std::vector<std::string>& a, \
    const std::vector<std::string>& b)
{
    std::vector<std::string> result;

    std::vector<std::string>::const_iterator it;
    for (it= a.begin(); it!=a.end(); ++it)
        if (std::find(b.begin(), b.end(), *it)==b.end())
            result.push_back(*it);

    return result;
}

///////////////////////////////////////////////////////////////////////////////
// Helper function update_ratemat (unfinished!!!)
template<class Type>
Eigen::SparseMatrix<Type> update_ratemat(
    const Eigen::SparseMatrix<Type>& ratemat, 
    const vector<Type>& state, 
    const vector<Type>& params, 
    const std::string& testwt_scale)
{
    Eigen::SparseMatrix<Type> result = ratemat;

    std::vector<int> inf_ind;
    int transm_ind;
    std::vector<int> transm_wt_ind;
    std::vector<int> foi_ind;
/*
    inf_ind = grep("I[a-z]", names(s)),
                            transm_ind = which(names(p) == "beta0"),
                            ## fragile! assumes same order as state
                            transm_wt_ind = grep("C[a-z]", names(p)),
                            foi_ind = c(which(rownames(M) == "S"),
                                        which(colnames(M) == "E"))
*/
    // We may need to check inf_ind, transm_ind, transm_wt_ind have valid values after
    // executing the following code block
    for (int i= 0; i< stateNames.size(); i++) 
        if (stateNames[i].at(0)=='I' && stateNames[i].at(1)>='a' && stateNames[i].at(1)<='z') 
            inf_ind.push_back(i);

    transm_ind = std::find(paramNames.begin(), paramNames.end(), "beta0") - paramNames.begin();
    std::cout << "transm_ind " << transm_ind << std::endl;
    
    for (int i= 0; i< paramNames.size(); i++)
        if (paramNames[i].at(0)=='C' && paramNames[i].at(1)>='a' && paramNames[i].at(1)<='z')
            transm_wt_ind.push_back(i);

    foi_ind.push_back(std::find(stateNames.begin(), stateNames.end(), "S") - stateNames.begin());
    foi_ind.push_back(std::find(stateNames.begin(), stateNames.end(), "E") - stateNames.begin());


    int ns = state.size();

    // update force of infection (foi) in rate matrix

    Type foi = 0;
    for (int i = 0; i < inf_ind.size(); i++) {
        // note TMB vectors use () not [] for indexing (but still zero-indexed)
        foi += state(inf_ind.at(i))*params(transm_wt_ind.at(i));
        std::cout << inf_ind.at(i) << ", " << transm_wt_ind.at(i) << std::endl;
    }
    std::cout << foi_ind.at(0) << ", " << foi_ind.at(1) << std::endl;
    foi *= params(transm_ind);

    result.coeffRef(foi_ind.at(0), foi_ind.at(1)) = foi;

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

/*
    result = 0; // just to make sure

    const Type* valPtr = mat.valuePtr();
    const int* rowIndexPtr = mat.innerIndexPtr();
    const int* outPtr = mat.outerIndexPtr();

    for (int j= 0; j<mat.outerSize(); j++)
        for (int i= outPtr[j]; i<outPtr[j+1]; i++) {
            int row = rowIndexPtr[i];
            result(row) = result(row) + valPtr[i];
        }
*/
    //std::cout << "ratemat = " << mat << std::endl;
    //std::cout << result << std::endl;

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

/*
    result = 0; // just to make sure

    const Type* valPtr = mat.valuePtr();
    const int* outPtr = mat.outerIndexPtr();

    for (int j= 0; j<mat.outerSize(); j++)
        for (int i= outPtr[j]; i<outPtr[j+1]; i++) {
            result(j) = result(j) + valPtr[i];
        }
*/
    //std::cout << "ratemat = " << mat << std::endl;
    //std::cout << result << std::endl;

    return result;
}

///////////////////////////////////////////////////////////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{
    struct timeval tv;

    gettimeofday(&tv, NULL);
    std::cout << "*** Begin of TMB_do_step() at " << tv.tv_sec*1000000+tv.tv_usec << std::endl;

    // hard-code the names of state and parameters in fixed order
    //std::vector<std::string> stateNames {"S", "E", "Ia", "Ip", "Im", "Is", "H", "H2", \
    //                            "ICUs", "ICUd", "D", "R", "X", "V"};
    // Perhaps we don't need parameter names, but I just keep it for possible future use
    //std::vector<std::string> paramNames {"beta0", "Ca", "Cp", "Cm", "Cs" \
    //                            "alpha", "sigma", "gamma_a", "gamma_m", "gamma_s", \
    //                            "gamma_p", "rho", "delta", "mu", "N" \
    //                            "E0", "nonhosp_mort", "iso_m", "iso_s", "phi1", \
    //                            "phi2", "psi1", "psi2", "psi3", "c_prop", \
    //                            "c_delay_mean", "c_delay_cv", "proc_disp", "zeta"};

    // Create a map (stateNames vector -> [0, 1, 2, ...])
    std::map<std::string, int> stateName2IndexMap;
    for (int i= 0; i< stateNames.size(); i++)
        stateName2IndexMap[stateNames[i]] = i; 

    // Create a map (paramNames vector -> [0, 1, 2, ...])
    std::map<std::string, int> paramName2IndexMap;
    for (int i= 0; i< paramNames.size(); i++)
        paramName2IndexMap[paramNames[i]] = i;


#ifdef TMB_VERBOSE
    std::map<std::string, int>::iterator it;
    for (it= stateName2IndexMap.begin(); it!=stateName2IndexMap.end(); ++it)
        std::cout << it->first << " => " << it->second  << std::endl;

    std::cout << "Ip" << " => " << stateName2IndexMap["Ip"] << std::endl;
#endif

    // Joint negative log-likelihood (stub)
    Type jnll= 0;

    // Get data and parameters from R
    DATA_VECTOR(state);
    DATA_SPARSE_MATRIX(ratemat);
    DATA_INTEGER(dt);			// integer->integer
    DATA_INTEGER(do_hazard);		// boolean->integer
    DATA_INTEGER(stoch_proc);		// boolean->integer
    DATA_INTEGER(do_exponential);	// boolean->integer
    DATA_STRING(testwt_scale);		// string->string

    PARAMETER_VECTOR(params);

#ifdef TMB_VERBOSE
    ratemat = ratemat * 100;
    std::cout << "state = " << state << std::endl;
    std::cout << "ratemat = " << ratemat << std::endl;
    std::cout << "dt = "<< dt << std::endl;
    std::cout << "do_hazard = " << do_hazard << std::endl;
    std::cout << "stoch_proc = " << stoch_proc << std::endl;
    std::cout << "do_exponential = " << do_exponential << std::endl;
    std::cout << "testwt_scale = "<< testwt_scale << std::endl;

    std::cout << "params = " << params << std::endl;
#endif

    // We've got everything we need, lets do the job ...

    std::vector<std::string> x_states {"X", "N", "P"};
    if (std::find(stateNames.begin(), stateNames.end(), "V") != stateNames.end())
        x_states.push_back("V");

    std::vector<std::string> p_states = exclude_states(stateNames, x_states);

    ratemat= update_ratemat(ratemat, state, params, testwt_scale);

    Eigen::SparseMatrix<Type> flows(ratemat.rows(), ratemat.cols()); // Just to get the right shape

    bool prpc_disp_exist = (std::find(paramNames.begin(), paramNames.end(), "proc_disp") != paramNames.end());

    if (!stoch_proc || (prpc_disp_exist && paramName2IndexMap["proc_disp"] < 0)) {
        if (!do_hazard) {
            flows= col_multiply(ratemat, state) * dt; // ## sweep(ratemat, state, MARGIN=1, FUN="*")*dt
            //std::cout << "ratemat = " << ratemat << std::endl;
            //std::cout << "state = " << state << std::endl;
            //std::cout << "dt = "<< dt << std::endl;
            //std::cout << "flows = " << flows << std::endl;
        } 
        else {
            vector<Type> S = rowSums(ratemat);
            vector<Type> E = exp(-S*dt);

            vector<Type> norm_sum = S;
            for (int i=0; i<norm_sum.size(); i++)
                if (norm_sum(i)!=0.0)
                    norm_sum(i) = state(i)/norm_sum(i);
            
            //std::cout << "S = " << S.size() << std::endl; 
            //std::cout << "S = " << S << std::endl;
            //std::cout << "E = " << -E << std::endl;

            E = 1 - E;
            flows = col_multiply(col_multiply(ratemat, norm_sum), E);
            //flows.diagonal().array() += 0; // this doesn't work before it requires the dianonal 
                                             // elements exist.
            for (int i=0; i<flows.rows(); i++)
                if (flows.coeff(i, i)!=0)
                    flows.coeffRef(i, i) = 0.0;
            //std::cout << "flows = " << flows << std::endl;
        }
    } else {
        //flows = 0; 
        //double t = rgamma(0.5, 1.0); //rgammawn(0.2, dt)
        // Below is the R code block
        /*
        flows <- ratemat ## copy structure
        flows[] <- 0
        for (i in seq(length(state))) {
            ## FIXME: allow Dirichlet-multinomial ?
            dW <- dt
            if (!is.na(proc_disp <- params[["proc_disp"]])) {
                dW <- pomp::rgammawn(sigma = proc_disp, dt = dt)
            }
            ## FIXME: need to adjust for non-conserving accumulators!
            flows[i, -i] <- pomp::reulermultinom(
                n = 1,
                size = state[[i]],
                rate = ratemat[i, -i],
                dt = dW
            )
        }
        */
    }

    vector<Type> outflow;
    if (!do_exponential) {
        Eigen::SparseMatrix<Type> flows1(flows.rows(), p_states.size());
        for (int i=0; i<p_states.size(); i++) {
            int colIndex = stateName2IndexMap[p_states[i]];
            flows1.col(i) = flows.col(colIndex);
        }
        outflow = rowSums(flows1);
        //outflow <- rowSums(flows[, p_states])
    } 
    else {
        // Below is the R code block
/*
        ## want to zero out outflows from S to non-S compartments
        ##  (but leave the inflows - thus we can't just zero these flows
        ##   out in the rate matrix!)
        S_pos <- grep("^S", rownames(ratemat), value = TRUE)
        notS_pos <- grep("^[^S]", colnames(ratemat), value = TRUE)
        notS_pos <- setdiff(notS_pos, x_states)
        outflow <- setNames(numeric(ncol(flows)), colnames(flows))
        ## only count outflows from S_pos to other S_pos (e.g. testing flows)
        outflow[S_pos] <- rowSums(flows[S_pos, S_pos, drop = FALSE])
        ## count flows from infected etc. to p_states (i.e. states that are *not* parallel accumulators)
        outflow[notS_pos] <- rowSums(flows[notS_pos, p_states])
*/
    }

    vector<Type> inflow = colSums(flows);

    //std::cout << "state = " << state << std::endl;
    state = state - outflow + inflow;
    //std::cout << "state = " << state << std::endl;

    // Below is the R code block
/*
    ## check conservation (*don't* check if we are doing an exponential sim, where we
    ##  allow infecteds to increase without depleting S ...)
    MP_badsum_action <- getOption("MP_badsum_action", "warning")
    MP_badsum_tol <- getOption("MP_badsum_tol", 1e-12)
    if (!do_exponential &&
        !(MP_badsum_action == "ignore") &&
        !stoch_proc) ## temporary: adjust reulermultinom to allow for x_states ...
        {
            calc_N <- sum(state[p_states])
            if (!isTRUE(all.equal(calc_N, sum(params[["N"]]), tolerance = MP_badsum_tol))) {
                msg <- sprintf("sum(states) != original N (delta=%1.2g)", sum(params[["N"]]) - calc_N)
                get(MP_badsum_action)(msg)
            }
        } ## not exponential run or stoch proc or ignore-sum

    return(state)
*/

    REPORT(state);

    gettimeofday(&tv, NULL);
    std::cout << "*** End of TMB_do_step() at   " << tv.tv_sec*1000000+tv.tv_usec << std::endl;

    std::cout << std::flush;

    return jnll;
}
   



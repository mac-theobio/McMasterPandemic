// update_ratemat2 is the function that applies our own "parser" to the formulas
// to calculate the non-zero elements in rate matrix.

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/time.h>
#include <TMB.hpp>

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
template<class Type>
bool has_testing(
    const vector<Type>* state,
    const vector<Type>* params= NULL,
    const Eigen::SparseMatrix<Type>* ratemat= NULL
)
{
    return false;
/*
    if (params!=NULL) {
        std::vector<std::string>::iteractor it = \
            std::find(paramNames.begin(), paramNames.end(), "testing_intensity");
        if (it!=paramNames.end())
            return (paramName2IndexMap["testing_intensity"]=0)
        else
            return true;
    }

    // if (ratemat) // in our simplified case, ratemat's row/col names are the same as state's names
     
        
has_testing <- function(state, params = NULL, ratemat = NULL) {
    if (!is.null(params)) {
        return("testing_intensity" %in% names(params) && params[["testing_intensity"]] != 0)
    }
    if (!is.null(ratemat)) {
        return(any(grepl("_t$", rownames(ratemat))))
    }
    return(any(grepl("_t$", names(state))))
*/
}

///////////////////////////////////////////////////////////////////////////////
template<class Type>
std::vector<int> pfun(
    std::string& from, 
    std:: string& to, 
    Eigen::SparseMatrix<Type>& mat, 
    bool value = FALSE, 
    bool recycle = FALSE) 
{
    std::vector<int> result;

    result.push_back(std::find(stateNames.begin(), stateNames.end(), from) - stateNames.begin());
    result.push_back(std::find(stateNames.begin(), stateNames.end(), to) - stateNames.begin());

    return result;

/*
    pfun_method <- getOption("macpan_pfun_method", "startsWith")
    find_pos_grep <- function(tag, x) {
        grep(sprintf("^%s(_|$)", tag), x, value = value)
    }
    find_pos_startsWith <- function(tag, x) {
        nt <- nchar(tag)
        r <- which(startsWith(x, tag)) ## subset quickly
        ## test remainder for _ or $
        r <- r[(substr(x[r], nt + 1, nt + 1) == "_") |
            nchar(x[r]) == nt]
        if (!value) {
            return(r)
        }
        return(x[r])
    }

    if (pfun_method == "both") {
        from_pos <- find_pos_grep(from, rownames(mat))
        to_pos <- find_pos_grep(to, colnames(mat))
        from_pos_sw <- find_pos_startsWith(from, rownames(mat))
        to_pos_sw <- find_pos_startsWith(to, colnames(mat))
        stopifnot(all(from_pos == from_pos_sw) && all(to_pos == to_pos_sw))
    } else {
        find_pos <- switch(pfun_method,
            startsWith = find_pos_startsWith,
            grep = find_pos_grep
        )
        from_pos <- find_pos(from, rownames(mat))
        to_pos <- find_pos(to, colnames(mat))
    }

    nf <- length(from_pos)
    nt <- length(to_pos)
    if (recycle && (nt == 1 || nf == 1)) {
        from_pos <- rep(from_pos, length.out = max(nt, nf))
        to_pos <- rep(to_pos, length.out = max(nt, nf))
    }
    if (!(length(to_pos) == length(from_pos) &&
        length(to_pos) > 0 && length(from_pos) > 0)) { ## must be positive
        stop(sprintf(
            "to_pos, from_pos don't match: from_pos=%s, to_pos=%s",
            paste(colnames(mat)[from_pos], collapse = ", "),
            paste(rownames(mat)[to_pos], collapse = ", ")
        ))
    }

    return(cbind(from_pos, to_pos))
*/
}

///////////////////////////////////////////////////////////////////////////////
template<class Type>
Type update_foi(
    vector<Type>& state, 
    vector<Type>& params, 
    vector<Type>& beta // *** This is a trouble-maker for R-C++ conversion as beta could be 
         // either a vector or a matrix. Here lets assume it is vector
) 
{
    assert(state.size()==beta.size());
    
    Type foi = 0;
    for (int i = 0; i < state.size(); i++) {
        foi += state(i)*beta(i);
    }

    return foi;

/*
    ## update infection rate
    if (is.matrix(beta)) {
        ## e.g. beta is a matrix when the transmission parameters are age-structured
        if (length(state) != ncol(beta)) stop("number of columns of beta must match length of state vector")
        foi <- beta %*% state
    } else {
        if (length(state) != length(beta)) {
            stop("length of state and beta are not the same")
        }
        foi <- sum(state * beta)
    }
    if (has_zeta(params)) {
        if (has_age(params)) stop("phenomenological heterogeneity (zeta != 0) untested with age-structed par
ams")
        ## suppose zeta is a vector zeta1, zeta2, zeta3, ...
        ##  we also need parameters   zeta_break1, zeta_break2 (0<zbi<1)
        ##  one *fewer* break parameter than zeta_i value
        ## if 0< S/N < zeta_break1   -> zeta1
        ##  zeta_break1 < S/N < zeta_break2 -> zeta2
        ## ...
        ##  zeta_breakx < S/N < 1  -> zetax
        Susc_frac <- 1 / params[["N"]] * sum(state[grep("^S_?", names(state))])
        if (any(grepl("zeta[0-9]", names(params)))) {
            zeta <- with(
                as.list(params),
                if (Susc_frac < zeta_break) zeta1 else zeta2
            )
        }
        ## alternately could just make it a vector
        ## ... but this messes with age-structured stuff
        ## if (length(zeta)>0) {
        ## ...
        ## }
        ## ???? het at pop level or sub-category level??
        foi <- foi * with(as.list(params), Susc_frac^zeta)
    }
    return(foi)
*/
}

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

    std::cout << "=========================================" << std::endl;
    std::cout << ratemat << std::endl;

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

    std::cout << "-------------------------------------" << std::endl;
    std::cout << result << std::endl;
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

    if (has_testing(&state, &params)) { // false in the simplied test case, so skip
        //aa <- c("wtsvec", "posvec", "testing_time")
        //saved_attrs <- setNames(lapply(aa, attr, x = ratemat), aa)
    }
    //## update testing flows. DO THIS FIRST, before updating foi: **** assignment via cbind() to Matrix objects loses attributes???
    if (has_testing(&state)) { // false in the simplied test case, so skip
        /*
        testing_time <- attr(ratemat, "testing_time") ## ugh  (see **** above)
        ## positions of untested, positive-waiting, negative-waiting compartments
        ## (flows from _n, _p to _t, or back to _u, are always at per capita rate omega, don't need
        ##  to be updated for changes in state)
        ## FIXME: backport to testify?
        u_pos <- grep("_u$", rownames(ratemat))
        p_pos <- grep("_p$", rownames(ratemat))
        n_pos <- grep("_n$", rownames(ratemat))
        ## original/unscaled prob of positive test by compartment, testing weights by compartment
        posvec <- attr(ratemat, "posvec")
        if (is.null(posvec)) stop("expected ratemat to have a posvec attribute")
        wtsvec <- attr(ratemat, "wtsvec")
        if (is.null(wtsvec)) stop("expected ratemat to have a wtsvec attribute")
        ## scaling ...
        testing_intensity <- params[["testing_intensity"]]
        testing_tau <- params[["testing_tau"]]
        N0 <- params[["N"]]
        W <- sum(wtsvec * state[u_pos])
        sc <- switch(testwt_scale,
            none = 1,
            N = N0 / W,
            sum_u = sum(state[u_pos]) / W,
            sum_smooth = {
                rho <- testing_intensity
                tau <- testing_tau
                tau * N0 / (tau * W + rho * N0)
                ## NOTE 'smoothing' doc has numerator rho*tau*N0,
                ## but testing intensity (rho) is included in ratemat
                ## calculation below ...
            }
        )
        ratemat[cbind(u_pos, n_pos)] <- testing_intensity * sc * wtsvec * (1 - posvec)
        ratemat[cbind(u_pos, p_pos)] <- testing_intensity * sc * wtsvec * posvec
        if (testing_time == "sample") {
            N_pos <- which(rownames(ratemat) == "N")
            P_pos <- which(rownames(ratemat) == "P")
            ratemat[cbind(u_pos, N_pos)] <- ratemat[cbind(u_pos, n_pos)]
            ratemat[cbind(u_pos, P_pos)] <- ratemat[cbind(u_pos, p_pos)]
        }
        */
    }

/*
    ## update FOIs
    ratemat[pfun("S", "E", ratemat)] <- update_foi(state, params, make_beta(state, params))

    ## update vaccine allocation rates
    if (has_vax(params)) {
        ratemat <- add_updated_vaxrate(state, params, ratemat)
    }

    ## ugh, restore attributes if necessary
    if (inherits(ratemat, "Matrix") && has_testing(params)) {
        for (a in aa) {
            attr(ratemat, a) <- saved_attrs[[a]]
        }
    }
    return(ratemat)
*/
    return result;

// Below is Ben's version of update_ratemat(...). Using it won't generate the same state as
// example("do_step")
/*
    Eigen::SparseMatrix<Type> result = ratemat;

    std::vector<int> inf_ind;
    int transm_ind;
    std::vector<int> transm_wt_ind;
    std::vector<int> foi_ind;

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
*/
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
    DATA_IVECTOR(from);
    DATA_IVECTOR(to);
    DATA_IVECTOR(count);
    DATA_IVECTOR(spi);
    DATA_IVECTOR(modifier);
    DATA_INTEGER(dt);			// integer->integer
    DATA_INTEGER(do_hazard);		// boolean->integer
    DATA_INTEGER(stoch_proc);		// boolean->integer
    DATA_INTEGER(do_exponential);	// boolean->integer
    DATA_STRING(testwt_scale);		// string->string

    PARAMETER_VECTOR(params);

    // Concatenate state and params
    vector<Type> sp(state.size()+params.size());
    sp << state, params;
    std::cout << "sp = " << sp << std::endl;

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

    ratemat= update_ratemat2(ratemat, sp, from, to, count, spi, modifier);

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

    vector<Type> inflow = colSums(flows);

    //std::cout << "state = " << state << std::endl;
    state = state - outflow + inflow;
    //std::cout << "state = " << state << std::endl;

    REPORT(state);

    gettimeofday(&tv, NULL);
    std::cout << "*** End of TMB_do_step() at   " << tv.tv_sec*1000000+tv.tv_usec << std::endl;

    std::cout << std::flush;

    return jnll;
}
   



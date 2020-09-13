## it takes existing ratemat from make_ratemat and transform it into the testify version rate mat

##' global variables for testify expansion
#' @export
non_expanded_states <- c("D","X")

##' @rdname non_expanded_states
##' @export
test_extensions <- c("u","p","n","t")

##' @rdname non_expanded_states
##' @export
test_accumulators <- c("N","P")

##' @rdname non_expanded_states
##' @export
asymp_cat <- c("S","E","Ia","Ip","R")

##' @rdname non_expanded_states
severe_cat <- c("Is","H","H2","ICUs","ICUd")

check_var_names <- function(var_names) {
    pat <- sprintf("^(%s)", paste(c("Im",asymp_cat,severe_cat), collapse="|"))
    extra_states <- grep(pat, var_names, pat, value=TRUE, invert=TRUE)
    if (length(extra_states)>0) {
        stop("states neither 'asymptomatic' nor 'severe' nor 'Im' :",paste(extra_states,collapse=", "))
    }
}

##' make weights vector for tests
##'
##' @param params parameter vector
##' @param var_names variables names, \emph{in matching order to state vector/rate matrix}
##' @examples
##' pp <- read_params("PHAC_testify.csv")
##' state1 <- state0 <- make_state(params=pp)   ## unexpanded
##' state1[] <- 1  ## occupy all states
##' state <- expand_stateval_testing(state0, params=pp)
##' vn <- setdiff(names(state0),non_expanded_states)
##' wtsvec <- make_test_wtsvec(pp, vn)
##' posvec <- make_test_posvec(pp, vn)
##' ## need to make_ratemate() with *unexpanded* state, then
##' ##  expand it
##' ratemat <- testify(make_ratemat(state1,pp), pp)
##' betavec <- make_betavec(state,pp)
##' tt2 <- ratemat
##' tt2[tt2>0] <- 1 ## make all edges == 1
##' show_ratemat(tt2)
##' show_ratemat(ratemat,aspect="fill")
##' if (require(igraph)) {
##'    g <- igraph::graph_from_adjacency_matrix(tt2)
##'    plot(g, layout=igraph::layout_nicely)
##' }
##' @export
make_test_wtsvec <- function(params,var_names=NULL) {
    if (!is.null(var_names)) {
        check_var_names(var_names)
    }
    ## test whether a specific set of W-parameters are the *only*
    ##  W-parameters in the parameter vector
    match_pars <- function(pars) {
        Wpars <- grep("^W",names(params),value=TRUE)
        return(identical(pars,sort(Wpars)))
    }
    if (match_pars("W_asymp")) {  ## W_asymp is the only weighting parameter
        ## one-parameter model: weights specified as (asymp=W_asymp, symp=1)
        symp_cat <- setdiff(var_names,asymp_cat)
        wts_vec <- rep(c(params[["W_asymp"]],1),
                       c(length(asymp_cat),length(symp_cat)))
        names(wts_vec) <- c(asymp_cat,symp_cat)
    } else if (match_pars(c("W_asymp","W_severe"))) {
        ## asymp= W_asymp (<1),  I_m=1, severe = W_severe (>1)
        mild_cat <- setdiff(var_names, c(asymp_cat, severe_cat))
        wts_vec <- rep(c(params[["W_asymp"]],1,params[["W_severe"]]),
                       c(length(asymp_cat),length(mild_cat),length(severe_cat)))
        names(wts_vec) <- c(asymp_cat,mild_cat,severe_cat)
    } else {
        ## general
        wts_vec <- params[grepl("^W",names(params))]
        names(wts_vec) <- gsub("^W","",names(wts_vec))
    }
    if (!all(sort(names(wts_vec))==sort(var_names))) {
        stop("weights vector names should match var names")
    }
    wts_vec <- wts_vec[var_names] ## reorder
    return(wts_vec)
}

##' make vector of test positivity
##'
##' @rdname make_test_wtsvec
##'
##' @param params parameter vector
##' @export
make_test_posvec <- function(params,var_names=NULL) {
    vec <- params[grepl("^P",names(params))]
    names(vec) <- gsub("^P_?","",names(vec))
    if (!is.null(var_names)) {
        if (!all(sort(names(vec))==sort(var_names))) {
            stop("vector names should match var names")
        }
        vec <- vec[var_names] ## reorder
    }
    return(vec)
}

##' expand states and values to include
##'
##' @param x state vector
##' @param method method for distributing values across new (expanded) states
##' @param add_accum add N and P (neg/pos test) accumulator categories?
##' @param params parameters
##' @examples
##' pp <- read_params("PHAC_testify.csv")
##' s <- make_state(params=pp)
##' expand_stateval_testing(s, params=pp)
##' expand_stateval_testing(s, method="untested")
##' @export

expand_stateval_testing <- function(x, method=c("eigvec","untested","spread"),
                            params=NULL,
                            add_accum=TRUE)
{
    method <- match.arg(method)
    
    newnames <- unlist(lapply(setdiff(names(x),non_expanded_states), paste, test_extensions, sep="_"))
    new_states <- rep(0,length(newnames))
    names(new_states) <- newnames
    ## FIXME: check on names, attributes, etc.
    if (method=="untested") {
        new_states[paste0(setdiff(names(x),non_expanded_states),"_u")] <- x[setdiff(names(x),non_expanded_states)]
    } else if (method=="spread") {
        ## slightly fragile: depends on ordering
        n_expand <- length(test_extensions)
        new_states <- rep(x, each=n_expand)/n_expand
    } else if (method=="eigvec") {
        S_tot <- x[["S"]]
        nonS_tot <- sum(x[names(x)!="S"])  ## or sum(x) - S_tot
        if (is.null(params)) stop("need params to expand via eigvec")
        ee <- rExp(params, testify=TRUE, return_val="eigenvector")
        ## watch out, ee doesn't include R ... match by name, not position
        is_S <- grep("^S", names(ee), value=TRUE)
        is_nonS <- grep("^S", names(ee), invert=TRUE, value=TRUE)
        new_states[is_S]  <- ee[is_S] *S_tot/sum(ee[is_S])
        new_states[is_nonS] <- ee[is_nonS]*nonS_tot/sum(ee[is_nonS])
    }
    if (add_accum) {
    ## FIXME, compare setNames = 0 code 
        new_states <- c(new_states, setNames(numeric(length(non_expanded_states)), non_expanded_states), c(N=0,P=0))
    }
    return(new_states)
}


##' expand rate matrix for testing status
##' @param ratemat original rate matrix
##' @param params parameters
##' @param testing_time "report" (N and P are counted at the time when individuals move from _n, _p to _u, _t compartments) or "sample" (N and P are counted when individuals move from _u to _n or _p)
##' @param debug what it sounds like
##' @examples
##' params <- read_params("PHAC_testify.csv")
##' state <- make_state(params[["N"]],E0=params[["E0"]])
##' M <- make_ratemat(state,params)
##' image(M)
##' @export
testify <- function(ratemat,params,debug=FALSE,
                    testing_time=NULL) {
    ## wtsvec is a named vector of testing weights which is the per capita rate at which ind in untested -> n/p
    ## truevec is a named vector of probability of true positive for (positive states), true negative for (negative states)
    ## Assuming false positive/negative is (1-trueprob)
    ## omega is waiting time for tests to be returned
    ## don't use match.args() because we may get handed NULL ...
    if (is.null(testing_time)) {
        testing_time <- "report"
    } else {
        if (!testing_time %in% c("report","sample")) {
            stop("unknown testing_time", testing_time)
        }
    }
    vn <- setdiff(rownames(ratemat), non_expanded_states)
    wtsvec <- make_test_wtsvec(params, var_names=vn)
    posvec <- make_test_posvec(params, var_names=vn)
    omega <- params[["omega"]]
    testing_intensity <- params[["testing_intensity"]]
    
    M <- ratemat  ## FIXME: why two names (M and ratemat)
    states <- rownames(ratemat)
    ## testifiable vars will get expanded
    expand_set <- setdiff(states, non_expanded_states)
    
    dummy_states <- setNames(numeric(length(expand_set)), expand_set)
    new_states <- names(expand_stateval_testing(dummy_states, method="untested"))
    
	ns <- length(new_states)
	new_M <- matrix(0,nrow=ns, ncol=ns
                 , dimnames=list(from=new_states,to=new_states)
                   )
   ## Between states
    for(i in rownames(ratemat)){
        sn <- function(state, compartment=i) paste0(compartment, "_", state)
        
        ## Between states: epidemiological transitions are the same as X -> Y
        for(j in colnames(ratemat)){
            if (debug) { 
      		print(cat(i,j,"\n"))
            }

            ## source expanded, destination expanded: (k_i -> k_j) at same rate as i -> j
            if((i %in% expand_set) && (j %in% expand_set)){
                for (k in test_extensions) {
                    new_M[sn(k,i),sn(k,j)] <- M[i,j]
                }
            }
            ## source expanded, destination not expanded: all k_i - > j at same rate as i -> j
            if((i %in% expand_set) && (j %in% non_expanded_states)){
                for (k in test_extensions) {
                    new_M[sn(k,i),j] <- M[i,j]
                }
            }
        }

        ## testing flows
   	if (i %in% expand_set) {
            new_M[sn("u"),sn("n")] <- testing_intensity*wtsvec[[i]]*(1-posvec[[i]])
            new_M[sn("u"),sn("p")] <- testing_intensity*wtsvec[[i]]*(posvec[[i]])
            new_M[sn("n"),sn("u")]  <- omega
            new_M[sn("p"),sn("t")]  <- omega
        } ## if i %in% expand_set
    }  ## loop over all states
    ## hospitalization special cases: everyone gets tested when they leave the Is compartment
    ##  for H, ICUs, or ICUd
    for (j in c("H","ICUs","ICUd")) {
        sn1 <- function(state, compartment) paste0(compartment, "_", state)
        new_M["Is_u",sn1("p",j)] <- new_M["Is_u",sn1("u",j)]*posvec[["Is"]]
        new_M["Is_u",sn1("n",j)] <- new_M["Is_u",sn1("u",j)]*(1-posvec[["Is"]])
        new_M["Is_u",sn1("u",j)] <- 0
    }
    if (inherits(M,"Matrix")) new_M <- Matrix(new_M)
    for (i in expand_set) {
        sn2 <- function(state, compartment=i) paste0(compartment, "_", state)
        if (testing_time=="report") {
            ## N, P are recorded at {n->u, p->t} transition (when tests are reported)
            new_M[sn2("n"),"N"] <- new_M[sn2("n"),sn2("u")] 
            new_M[sn2("p"),"P"] <- new_M[sn2("p"),sn2("t")] 
        } else {
            ## N, P are recorded at {u->n, u->p} transition (when samples are taken)
            new_M[sn2("u"),"N"] <- new_M[sn2("u"),sn2("n")]
            new_M[sn2("u"),"P"] <- new_M[sn2("u"),sn2("p")] 
        }
    }

    attr(new_M,"wtsvec") <- wtsvec
    attr(new_M,"posvec") <- posvec
    attr(new_M,"testing_time") <- testing_time
    return(new_M)
}

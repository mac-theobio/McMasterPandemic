## it takes existing ratemat from make_ratemat and transform it into the testify version rate mat

##' make weights vector for tests
##'
##' @param params parameter vector
##' @param var_names variables names, \emph{in matching order to state vector/rate matrix}
##' @examples
##' pp <- read_params("PHAC_testify.csv")
##' state1 <- state0 <- make_state(params=pp)   ## unexpanded
##' state1[] <- 1  ## occupy all states
##' state <- expand_stateval(state0)
##' wtsvec <- make_test_wtsvec(pp, var_names=setdiff(names(state0),c("D","R","X")))
##' posvec <- make_test_posvec(pp)
##' ## need to make_ratemate() with *unexpanded* state, then
##' ##  expand it
##' ratemat <- testify(make_ratemat(state1,pp), pp)
##' betavec <- make_betavec(state,pp)
##' tt2 <- ratemat
##' tt2[tt2>0] <- 1 ## make all edges == 1
##' if (require(Matrix)) {
##'     Matrix::image(Matrix(tt2))
##' }
##' heatmap(ratemat, Rowv=NA, Colv=NA, scale="none")
##' if (require(igraph)) {
##'    g <- igraph::graph_from_adjacency_matrix(tt2)
##'    plot(g, layout=igraph::layout_nicely)
##' }
##' @export
make_test_wtsvec <- function(params,var_names=NULL) {
    ## FIXME: do normalization here??
    ## only need var_names if we're going to set by asymp/symp
    if (identical("W_asymp",grep("^W",names(params),value=TRUE))) {
        W_asymp <- params[["W_asymp"]]
        ## one-parameter model: weights specified as (W_asymp, 1-W_asymp)
        asymp_cat <- c("S","E","Ia")
        symp_cat <- setdiff(var_names,asymp_cat)
        wts_vec <- setNames(rep(c(W_asymp,1-W_asymp),
                                c(length(asymp_cat),length(symp_cat))),
                            c(asymp_cat,symp_cat))
    } else {
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
make_test_posvec <- function(params) {
    pos_vec <- params[grepl("P",names(params))]
    return(pos_vec)
}

# ## use inside testify to expand states
# ## DRY: use this inside expand_stateval??
# expand_states <- function(expandable,nonexpandable) {
#     new_states <- c(paste0(expandable,c("_u","_p","_n","_t")),nonexpandable, "N", "P")
# }

##' expand states and values to include
##'
##' @param x state vector
##' @param method method for distributing values across new (expanded) states
##' @param add_accum add N and P (neg/pos test) accumulator categories?
##' @param non_expanded states that should \emph{not} be expanded into testing categorieso
##' @export
expand_stateval <- function(x, method=c("untested","spread"),
                            add_accum=TRUE,
                            non_expanded = c("R","D","X")){
    method <- match.arg(method)
    extensions <- c("u","p","n","t")
    newnames <- unlist(lapply(setdiff(names(x),non_expanded), paste, extensions, sep="_"))
    new_states <- rep(0,length(newnames))
    names(new_states) <- newnames
    if (method=="untested") {
        new_states[paste0(setdiff(names(x),non_expanded),"_u")] <- x[setdiff(names(x),non_expanded)]
    } else {
        ## slightly fragile: depends on ordering
        n_expand <- length(extensions)
        new_states <- rep(x, each=n_expand)/n_expand
    }
    if (add_accum) {
    ## FIXME, compare setNames = 0 code 
        new_states <- c(new_states, setNames(numeric(length(non_expanded)), non_expanded), c(N=0,P=0))
    }
    return(new_states)
}


##' expand rate matrix for testing status
##' @param ratemat original rate matrix
##' @param params parameters
##' @param non_expand_set states \emph{not} to expand
##' @param debug what it sounds like
##' @export
testify <- function(ratemat,params,debug=FALSE,
                    non_expand_set=c("D","R","X")) {
    ## wtsvec is a named vector of testing weights which is the per capita rate at which ind in untested -> n/p
    ## truevec is a named vector of probability of true positive for (positive states), true negative for (negative states)
    ## Assuming false positive/negative is (1-trueprob)
    ## omega is waiting time for tests to be returned
    wtsvec <- make_test_wtsvec(params, var_names=setdiff(rownames(ratemat),
                                                         non_expand_set))
    posvec <- make_test_posvec(params)
    omega <- params[["omega"]]
    testing_intensity <- params[["testing_intensity"]]
    
    M <- ratemat  ## FIXME: why two names?
    states <- rownames(ratemat)
    ## testifiable vars will get expanded
    ## We need to change expand_set here!! 

   ## Following page 8 L137 in McMasterReport2020-07-06.pdf 
    ## expand_set <- c("S","E","Ia","Ip","Im","Is","H","H2","hosp","ICUs","ICUd","D","R","X")
    expand_set <- setdiff(states, non_expand_set)
    ## If these are not expandable, we need to get rid of WD, WR, WX, PD, PR, PX
    ## non_expand_set <- c("P","N") ## adding P and N for accumulation
    dummy_states <- setNames(numeric(length(expand_set)), expand_set)
    new_states <- names(expand_stateval(dummy_states))
    
	ns <- length(new_states)
	new_M <- matrix(0,nrow=ns, ncol=ns
                 , dimnames=list(from=new_states,to=new_states)
                   )
   ## Between states
    for(i in rownames(ratemat)){
        ## Between states: horizontal arrows moves the same rate
        for(j in colnames(ratemat)){
            if (debug) { 
      		print(cat(i,j,"\n"))
            }
            if((i %in% expand_set) && (j %in% expand_set)){ 
           	new_M[paste0(i,"_u"),paste0(j,"_u")] <- M[i,j]
                new_M[paste0(i,"_t"),paste0(j,"_t")] <- M[i,j]
                new_M[paste0(i,"_p"),paste0(j,"_p")] <- M[i,j] 
                new_M[paste0(i,"_n"),paste0(j,"_n")] <- M[i,j]
            }
            if((i %in% expand_set) && (j %in% non_expand_set)){
                new_M[paste0(i,"_u"),j] <- M[i,j]
                new_M[paste0(i,"_t"),j] <- M[i,j]
                new_M[paste0(i,"_p"),j] <- M[i,j]
                new_M[paste0(i,"_n"),j] <- M[i,j]
            }
        }

        pn <- function(par, compartment=i) paste0(par, compartment)
        sn <- function(state, compartment=i) paste0(compartment, "_", state)
        
   	if (i %in% expand_set){
            new_M[sn("u"),sn("p")] <- testing_intensity*wtsvec[i]*(posvec[pn("P")])
            new_M[sn("u"),sn("n")] <- testing_intensity*wtsvec[i]*(1-posvec[pn("P")])
            new_M[sn("n"),sn("u")] <- new_M[sn("n"),"N"] <- omega
            new_M[sn("p"),sn("t")] <- new_M[sn("p"),"P"] <- omega
   	}
    }
    attr(new_M,"wtsvec") <- wtsvec
    attr(new_M,"posvec") <- posvec
    return(new_M)
}

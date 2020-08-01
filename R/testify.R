## it takes existing ratemat from make_ratemat and transform it into the testify version rate mat



##' make weights vector for tests
##'
##' @param params parameter vector
##' 
##' @examples
##' pp <- read_params("PHAC_testify.csv")
##' state1 <- state0 <- make_state(params=pp)   ## unexpanded
##' state1[] <- 1  ## occupy all states
##' state <- expand_stateval(state0)
##' wtsvec <- make_test_wtsvec(pp)
##' posvec <- make_test_posvec(pp)
##' ## need to make_ratemate() with *unexpanded* state, then
##' ##  expand it
##' ratemat <- testify(make_ratemat(state1,pp), pp)
##' betavec <- make_betavec(state,pp)
##' tt2 <- ratemat
##' tt2[tt2>0] <- 1 ## make all edges == 1
##' Matrix::image(Matrix(tt2))
##' heatmap(ratemat, Rowv=NA, Colv=NA, scale="none")
##' if (require(igraph)) {
##'    g <- igraph::graph_from_adjacency_matrix(tt2)
##'    plot(g, layout=igraph::layout_nicely)
##' }
##' @export
make_test_wtsvec <- function(params) {
    ## FIXME: why is this hard coded rather than pasted?
    wtscats <- c("WS","WE","WIa","WIp","WIm","WIs","WH","WH2","Whosp","WICUs","WICUd","WD","WR","WX")
    wts_vec <- with(as.list(params),setNames(c(WS,WE,WIa,WIp,WIm,WIs,WH,WH2,Whosp,WICUs,WICUd,WD,WR,WX),wtscats))
    return(wts_vec)
}

update_testflows <- function(state, params, wts_vec) {
}

##' make vector of test positivity
##'
##' @rdname make_test_wtvec
##'
##' @param params parameter vector
##' @export
make_test_posvec <- function(params) {
    ## FIXME: why is this hard coded rather than pasted?
    poscats <- c("PS","PE","PIa","PIp","PIm","PIs","PH","PH2","Phosp","PICUs","PICUd","PD","PR","PX")
    pos_vec <- with(as.list(params),setNames(c(PS,PE,PIa,PIp,PIm,PIs,PH,PH2,Phosp,PICUs,PICUd,PD,PR,PX),poscats))
    return(pos_vec)
}

## use inside testify to expand states
## DRY: use this inside expand_stateval??
expand_states <- function(x) {
    new_states <- c(paste0(x,c("_u","_p","_n","_t")), "N", "P")
}

##' expand states and values to include
##'
##' @param x state vector
##' @param method method for distributing values across new (expanded) states
##' @param add_accum add N and P (neg/pos test) accumulator categories?
##' 
##' @export
expand_stateval <- function(x, method=c("untested","spread"),
                            add_accum=TRUE) {
    method <- match.arg(method)
    extensions <- c("u","p","n","t")
    newnames <- unlist(lapply(names(x), paste, extensions, sep="_"))
    new_states <- rep(0,length(newnames))
    names(new_states) <- newnames
    if (method=="untested") {
        new_states[paste0(names(x),"_u")] <- x
    } else {
        ## slightly fragile: depends on ordering
        n_expand <- length(extensions)
        new_states <- rep(x, each=n_expand)/n_expand
    }
    if (add_accum) {
        new_states <- c(new_states, c(N=0,P=0))
    }
    return(new_states)
}


##' expand rate matrix for testing status
##' @param ratemat original rate matrix
##' @param params parameters
##' @param debug what it sounds like
##' @export
testify <- function(ratemat,params,debug=FALSE){
    ## wtsvec is a named vector of testing weights which is the per capita rate at which ind in untested -> n/p
    ## truevec is a named vector of probability of true positive for (positive states), true negative for (negative states)
    ## Assuming false positive/negative is (1-trueprob)
    ## omega is waiting time for tests to be returned
    wtsvec <- make_test_wtsvec(params)
    posvec <- make_test_posvec(params)
    omega <- params[["omega"]]
	
    M <- ratemat  ## FIXME: why two names?
    states <- rownames(ratemat)
    ## testifiable vars will get expanded

   ## Following page 8 L137 in McMasterReport2020-07-06.pdf 
    ## expand_set <- c("S","E","Ia","Ip","Im","Is","H","H2","hosp","ICUs","ICUd","D","R","X")
    ## non_expand_set <- c("P","N") ## adding P and N for accumulation
    dummy_states <- setNames(numeric(length(states)), states)
    new_states <- names(expand_stateval(dummy_states))
    
   ns <- length(new_states)
   new_M <- matrix(0,nrow=ns, ncol=ns
                 , dimnames=list(from=new_states,to=new_states)
                   )

    expand_set <- setdiff(states, c("D","R","X"))
	
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
        }
	new_M[paste0(i,"_u"),paste0(i,"_p")] <- wtsvec[paste0("W",i)]*(posvec[paste0("P",i)])
	new_M[paste0(i,"_u"),paste0(i,"_n")] <- wtsvec[paste0("W",i)]*(1-posvec[paste0("P",i)])
	new_M[paste0(i,"_n"),paste0(i,"_u")] <- new_M[paste0(i,"_n"),"N"] <- omega
	new_M[paste0(i,"_p"),paste0(i,"_t")] <- new_M[paste0(i,"_p"),"P"] <- omega
	
    }
    return(new_M)
}

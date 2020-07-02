## it takes existing ratemat from make_ratemat and transform it into the testify version rate mat

##' @examples
##' library(McMasterPandemic)
##' p <- read_params("ICU1.csv")
##' state <- make_state(params=p)
##' ratemat <- make_ratemat(state,p)
##' tt <- testify(ratemat, debug=TRUE)
##' tt2 <- tt
##' tt2[tt2>0] <- 1 ## make all edges == 1 
##' if (require(igraph)) {
##'    g <- igraph::graph_from_adjacency_matrix(tt2)
##'    plot(g, layout=igraph::layout_nicely)
##' }
testify <- function(ratemat,test_S=0.5, test_E=0.5, test_Ia=0.5, test_Ip=0.5, test_Im = 0.5, test_Is=0.5,
                    omega=0.5, debug=FALSE){
	 ## omegaT = rate of untested to testing?? X_u -> X_t
	 ## omega = rate of testing to untested, also rate of moving to N; X_t -> X_u and X_t -> N
	 ## omega = rate of testing to positive, also rate of moving to P; X_t -> X_p and X_t -> P
    ## categorizing states 
    negative <- c("S","E")
    ## BMB: I wouldn't bother including the non-testifiable vars in 'positive'???
    positive <- c("Ia","Ip","Im","Is")
    testifiable_var <- c(negative,positive)
	
    M <- ratemat
    states <- rownames(ratemat)
	
    ## testifiable vars will get expanded
	
    expand_states <- function(x) {
        new_states <- x
        if(x %in% testifiable_var){
            new_states <- paste0(x,c("_u","_t"))
            if(x %in% positive){
                new_states <- c(new_states,paste0(x,"_p"))
            }
        }
        return(new_states)
    }

    ## BMB: can we skip non-testifiable states entirely?
    ## MLi: We need non-testifiable states for j
    testified_states <- unlist(lapply(states,expand_states))
    testified_states <- c(testified_states,"H","hosp","X","ICUs","ICUd","H2","D","P", "N") ## accumulate P and N and adding back all the non-testify states
	
    ns <- length(testified_states)
    new_M <- matrix(0,nrow=ns, ncol=ns
                  , dimnames=list(from=testified_states,to=testified_states))
	
    for(i in rownames(ratemat)){
        for(j in colnames(ratemat)){
            if (debug) print(cat(i,j,"\n"))
            if( (i!=j) && (M[i,j] > 0)){  ## avoid R -> X and D -> X
                ## Both vars are not testifiable var
                if(!(i %in% testifiable_var) & !(j %in% testifiable_var)){
                    new_M[i,j] <- M[i,j]
                }
                ## testifiable to not testifiable
                if((i %in% testifiable_var) & !(j %in% testifiable_var)){
                    if(j != "R") { 
                        new_M[paste0(i,"_u"),j] <- new_M[paste0(i,"_u"),"P"] <- M[i,j]   
                        new_M[paste0(i,"_t"),j] <- new_M[paste0(i,"_t"),"P"] <- M[i,j]
                        new_M[paste0(i,"_p"),j] <- M[i,j]
                    }
                    if(j == "R"){ new_M[paste0(i,"_u"),"P"] <- new_M[paste0(i,"_t"),"P"] <- 0} ## dont need?
                }
                ## Both testifiable (both negative or positive will have the same 2 transitions)
                if(((i %in% negative) & (j %in% negative))
                		| ((i %in% positive) & (j %in% positive))){
                    new_M[paste0(i,"_u"),paste0(j,"_u")] <- M[i,j]
                    new_M[paste0(i,"_t"),paste0(j,"_t")] <- M[i,j]
                }
                    ## BMB: this is different for neg and pos states 
                    ## MLi: search "N" and "P"
                    ## neg states: add a flow to N
                    ## pos states: add a flow to P
                    ## positive to positive gets an additional transition
                    if((i %in% positive) & (j %in% positive)){
                        new_M[paste0(i,"_p"),paste0(j,"_p")] <- M[i,j]
                    }
                    ## Po
                    if((i %in% negative) & (j %in% positive)){ ## cross flow from E_t -> Ia_u or Ip_u
                    	   new_M[paste0(i,"_u"),paste0(j,"_u")] <- M[i,j]
                    	   new_M[paste0(i,"_t"),paste0(j,"_u")] <- M[i,j]


                    }
                }
            }
        }
        ## within state
        for(i in rownames(ratemat)){
            if(i %in% testifiable_var){
                if(i %in% negative){
                    new_M[paste0(i,"_t"),paste0(i,"_u")] <- new_M[paste0(i,"_t"),"N"] <- omega
                }
                if(i %in% positive){
                    new_M[paste0(i,"_t"),paste0(i,"_p")] <- new_M[paste0(i,"_t"),"P"] <- omega
                }
            }
        }
    ## manually inputing different test prop for each testifiable state
    new_M["S_u","S_t"] <- test_S
    new_M["E_u","E_t"] <- test_E
    new_M["Ia_u","Ia_t"] <- test_Ia
    new_M["Ip_u","Ip_t"] <- test_Ip
    new_M["Im_u","Im_t"] <- test_Im
    new_M["Is_u","Is_t"] <- test_Is
    return(new_M)
}

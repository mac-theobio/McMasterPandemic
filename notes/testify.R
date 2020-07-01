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
testify <- function(ratemat,eta=0.5,lambda=0.5, omegaT=0.5,
                    omegaN=0.5,
                    omegaP=0.5, debug=FALSE){
    ## BMB: describe lambda, omegaN, omegaP, omegaT?
    ## categorizing states 
    negative <- c("S","E")
    ## BMB: I wouldn't bother including the non-testifiable vars in 'positive'???
    positive <- c("Ia","Ip","Im","Is","H","hosp","X","ICUs","ICUd","H2","D")
    testifiable_var <- c("S","E","Ia","Ip","Im","Is")
	
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
    testified_states <- unlist(lapply(states,expand_states))
    testified_states <- c(testified_states, "P", "N") ## accumulate P and N
	
    ns <- length(testified_states)
    new_M <- matrix(0,nrow=ns, ncol=ns
                  , dimnames=list(from=testified_states,to=testified_states))
	
    for(i in rownames(ratemat)){
        for(j in colnames(ratemat)){
            if (debug) cat(i,j,"\n")
            if( (i!=j) && (M[i,j] > 0)){  ## avoid R -> X and D -> X
                ## Both vars are not testifiable var
                if(!(i %in% testifiable_var) & !(j %in% testifiable_var)){
                    new_M[i,j] <- M[i,j]
                }
                ## testifiable to not testifiable
                if((i %in% testifiable_var) & !(j %in% testifiable_var)){
                    if(j != "R") { 
                        new_M[paste0(i,"_u"),j] <- new_M[paste0(i,"_u"),"P"] <- (1-eta)*M[i,j]   ## eta and lambda are temporary the same across states, we can be more flexible later
                        new_M[paste0(i,"_t"),j] <- new_M[paste0(i,"_t"),"P"] <-eta*(1-lambda)*M[i,j]
                        new_M[paste0(i,"_p"),j] <- eta*(lambda)*M[i,j]
                    }
                    if(j == "R"){ new_M[paste0(i,"_u"),"P"] <- new_M[paste0(i,"_t"),"P"] <- 0}
                }
                ## Both testifiable (both negative and positive will have the same 3 transitions)
                if((i %in% testifiable_var) & (j %in% testifiable_var)){
                    new_M[paste0(i,"_u"),paste0(j,"_u")] <- (1-eta)*M[i,j]
                    new_M[paste0(i,"_t"),paste0(j,"_t")] <- (1-eta)*(1-lambda)*M[i,j]
                    ## BMB: this is different for neg and pos states
                    ## neg states: add a flow to N
                    ## pos states: add a flow to P
                    ## new_M[paste0(i,"_t"),paste0(j,"_u")] <- XX # tested but moved before result? Cross flow
                    ## positive to positive gets an additional transition
                    if((i %in% positive) & (j %in% positive)){
                        new_M[paste0(i,"_p"),paste0(j,"_p")] <- (1-eta)*(lambda)*M[i,j]
                    }
                }
            }
        }
        ## within state
        for(i in rownames(ratemat)){
            if(i %in% testifiable_var){
                new_M[paste0(i,"_u"),paste0(i,"_t")] <- omegaT
                if(i %in% negative){
                    new_M[paste0(i,"_t"),paste0(i,"_u")] <- new_M[paste0(i,"_t"),"N"] <- omegaN
                }
                if(i %in% positive){
                    new_M[paste0(i,"_t"),paste0(i,"_p")] <- new_M[paste0(i,"_t"),"P"] <- omegaP
                }
            }
        }
    }
    return(new_M)
}

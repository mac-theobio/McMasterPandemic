## it takes existing ratemat from make_ratemat and transform it into the testify version rate mat

##' @examples
##' library(McMasterPandemic)
##' p <- read_params("ICU1.csv")
##' state <- make_state(params=p)
##' wts <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,0.91,0.92,0.93,0.94)  ## This should be in ICU1.csv with all the other parameters
##' names(wts) <- c("WS","WE","WIa","WIp","WIm","WIs","WH","WH2","Whosp","WICUs","WICUd","WD","WR","WX")
##' posprob <- rep(c(0.05,0.9),c(2,12))
##' names(posprob) <- c("PS","PE","PIa","PIp","PIm","PIs","PH","PH2","Phosp","PICUs","PICUd","PD","PR","PX")
##' ratemat <- make_ratemat(state,p)
##' tt <- testify(ratemat,wts,posprob,omega=0.1,debug=TRUE)
##' tt2 <- tt
##' tt2[tt2>0] <- 1 ## make all edges == 1 
##' if (require(igraph)) {
##'    g <- igraph::graph_from_adjacency_matrix(tt2)
##'    plot(g, layout=igraph::layout_nicely)
##' }
##' 
##' 
##' 
##' 

make_test_wtsvec <- function(params) {
    wtscats <- c("WS","WE","WIa","WIp","WIm","WIs","WH","WH2","Whosp","WICUs","WICUd","WD","WR","WX")
    wts_vec <- with(as.list(params),setNames(c(WS,WE,WIa,WIp,WIm,WIs,WH,WH2,Whosp,WICUs,WICUd,WD,WR,WX),wtscats))
    return(wts_vec)
}

make_test_posvec <- function(params) {
    poscats <- c("PS","PE","PIa","PIp","PIm","PIs","PH","PH2","Phosp","PICUs","PICUd","PD","PR","PX")
    pos_vec <- with(as.list(params),setNames(c(PS,PE,PIa,PIp,PIm,PIs,PH,PH2,Phosp,PICUs,PICUd,PD,PR,PX),truecats))
    return(pos_vec)
}


testify <- function(ratemat,wtsvec,posvec,omega, debug=FALSE){
	# wtsvec is a named vector of testing weights which is the per capita rate at which ind in untested -> n/p
	# truevec is a named vector of probability of true positive for (positive states), true negative for (negative states)
	# Assuming false positive/negative is (1-trueprob)
	# omega is waiting time for tests to be returned
	M <- ratemat
   states <- rownames(ratemat)
   ## testifiable vars will get expanded
	
   expand_states <- function(x) {
   	new_states <- paste0(x,c("_u","_p","_n","_t"))
   }
   ## Following page 8 L137 in McMasterReport2020-07-06.pdf 
   expand_set <- c("S","E","Ia","Ip","Im","Is","H","H2","hosp","ICUs","ICUd","D","R","X")
   non_expand_set <- c("P","N") ## adding P and N for accumulation
   new_states <- c(unlist(lapply(expand_set,expand_states)), non_expand_set)

    
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
         if((i %in% expand_set) & (j %in% expand_set)){ 
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
  
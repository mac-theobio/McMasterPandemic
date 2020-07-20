## it takes existing ratemat from make_ratemat and transform it into the testify version rate mat

##' @examples
##' library(McMasterPandemic)
##' p <- read_params("ICU1.csv")
##' state <- make_state(params=p)
##' params <- c(0.1,0.2,0.3,0.4,0.5,0.6)  ## This should be in ICU1.csv with all the other parameters
##' names(params) <- c("WS","WE","WIa","WIp","WIm","WIs")
##' ratemat <- make_ratemat(state,p)
##' tt <- testify(ratemat, params,omega=0.1,debug=TRUE)
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

make_testvec <- function(params) {
    Testcats <- c("WS","WE","WIa","WIp","WIm","WIs")
    ## NB meaning of iso_* has switched from Stanford model
    test_vec0 <- with(as.list(params),setNames(c(Wsusp,We,Wa,Wp,Wm,Ws),Testcats))
    return(test_vec)
}


testify <- function(ratemat,testvec,omega, debug=FALSE){
	# testvec is a named vector of testing weights which is the per capita rate at which ind in untested -> n/p
	# omega is waiting time for tests to be returned
	M <- ratemat
   states <- rownames(ratemat)
   ## testifiable vars will get expanded
	
   expand_states <- function(x) {
   	new_states <- paste0(x,c("_u","_p","_n","_t"))
   }
   ## Following page 8 L137 in McMasterReport2020-07-06.pdf 
   expand_set <- c("S","E","Ia","Ip","Im","Is")
   non_expand_set <- c("H","H2","hosp","ICUs","ICUd","D","R","X","P","N") ## adding P and N for accumulation
   new_states <- c(unlist(lapply(expand_set,expand_states)), non_expand_set)

    
   ns <- length(new_states)
   new_M <- matrix(0,nrow=ns, ncol=ns
		, dimnames=list(from=new_states,to=new_states)
	 	)
	
   ## Between states
	for(i in rownames(ratemat)){
		for(j in colnames(ratemat)){
      	if (debug) print(cat(i,j,"\n"))
        	## Both testifiable
         	if((i %in% expand_set) & (j %in% expand_set)){ ## horizontal arrows moves the same rate
           		new_M[paste0(i,"_u"),paste0(j,"_u")] <- M[i,j]
            	new_M[paste0(i,"_t"),paste0(j,"_t")] <- M[i,j]
            	new_M[paste0(i,"_p"),paste0(j,"_p")] <- M[i,j] ## except for E -> I
            	if(i == "E"){
            	new_M[paste0(i,"_p"),paste0(j,"_p")] <- 0
            }
        		if((i == "S")&(j == "E")){
            	new_M[paste0(i,"_n"),paste0(j,"_n")] <- M[i,j] 
            	new_M[paste0(i,"_p"),paste0(j,"_p")] <- 0
        		}
            }
        	   ## Both vars are not testifiable var
            if(!(i %in% expand_set) & !(j %in% expand_set)){
                    new_M[i,j] <- M[i,j]
            }
            ## testifiable to not testifiable
            if((i %in% expand_set) & !(j %in% expand_set)){
            	if(j != "R") { 
               	new_M[paste0(i,"_u"),j] <- new_M[paste0(i,"_u"),"P"] <- M[i,j]   
                  new_M[paste0(i,"_t"),j] <- new_M[paste0(i,"_t"),"P"] <- M[i,j]
                  new_M[paste0(i,"_p"),j] <- M[i,j]
            	}
          	}
		}
	}
               
        ## within state
	for(i in expand_set){
   	if(i %in% c("S","E")){
      	new_M[paste0(i,"_u"),paste0(i,"_n")] <- testvec[paste0("W",i)]
         new_M[paste0(i,"_n"),paste0(i,"_u")] <- new_M[paste0(i,"_n"),"N"] <- omega
     	}
      if(i %in% c("Ia","Ip","Im","Is")){
      	new_M[paste0(i,"_u"),paste0(i,"_p")] <- testvec[paste0("W",i)]
         new_M[paste0(i,"_p"),paste0(i,"_t")] <- new_M[paste0(i,"_p"),"P"] <- omega
      }
	}
   return(new_M)
}
  

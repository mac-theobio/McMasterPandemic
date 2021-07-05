library(mvbutils)
library(igraph)
library(McMasterPandemic)
library(visNetwork)
pos <- which(search()=="package:McMasterPandemic")
## FIXME: include predict methods?
funs <- grep("mle_fun|forecast|calibrate|run_|do_step",ls(pos=pos),value=TRUE)
funs <- c(funs,"predict.fit_pansim")
ff <- foodweb(where=pos,
              funs=funs,
              ## rprune="forecast|calibrate|run_",
              plotting=FALSE, ancestors=FALSE,descendents=FALSE)
## HACK: foodweb doesn't recognize do.call() or ... ?
M <- ff$funmat
M["predict.fit_pansim",c("forecast_sim","forecast_ensemble")] <- 1
## HACK: set up another function parallel to run_sim_break
## M <- rbind(M,0)
## M <- cbind(M,0)
## rownames(M)[nrow(M)] <- "OTHER_TIMEFUN"
## colnames(M)[ncol(M)] <- "OTHER_TIMEFUN"
## M[,"OTHER_TIMEFUN"] <- M[,"run_sim_break"]
## M["OTHER_TIMEFUN",] <- M["run_sim_
## HACK: calibrate effectively calls the other run_sim
run_sim_funs <- setdiff(grep("run_sim_",rownames(M),value=TRUE), "run_sim_range")
M["forecast_sim",run_sim_funs] <- 1
M[run_sim_funs,"run_sim"] <- 1
g <- igraph::graph_from_adjacency_matrix(M)
plot(g,layout=layout.graphopt) ## fruchterman.reingold
e <- which(M==1, arr.ind=TRUE)
e[] <- rownames(M)[e]
e <- setNames(as.data.frame(e), c("from", "to"))
visNetwork(nodes=data.frame(id=rownames(M),
                            label=rownames(M)),
           edges=e) %>%
  visEdges(arrows = "to") %>%
  visHierarchicalLayout(sortMethod = "directed")

##https://datastorm-open.github.io/visNetwork/layout.html

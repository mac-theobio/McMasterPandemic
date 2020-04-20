## FIXME: clean up/minimize

##' display flow chart on current graphics device
##' @importFrom igraph graph_from_adjacency_matrix layout_as_tree
##' @importFrom Matrix Matrix
##' @importFrom graphics image
##' @param params parameters (flowchart is determined by non-zero flows)
##' @export
make_flowchart <- function(params=read_params("ICU1.csv")) {
    ## require(Matrix)
    state <- make_state(N=1e6,E0=1)
    state[] <- 1  ## all population states occupied
    M <- make_ratemat(state,params,do_ICU=TRUE)
    M2 <- M
    M2[M2>0] <- 1 ## make all edges == 1 (otherwise edges are so narrow they disappear; could scale up or plot as flows instead if we wanted to show strength?
    ## image(Matrix(M2))
    g <- igraph::graph_from_adjacency_matrix(M2)
    plot(g, layout=igraph::layout_as_tree)
}

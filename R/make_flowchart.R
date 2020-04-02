## FIXME: clean up/minimize

##' @importFrom igraph graph_from_adjacency_matrix
##' @importFrom Matrix Matrix
make_flowchart <- function() {
    params <- read_params(system.file("params","ICU1.csv",package="McMasterPandemic"))
    sdate <- "10-Feb-2020" ## arbitrary!
    state <- c(S=params[["N"]],E=params[["E0"]],Ia=0,Ip=0,Im=0,Is=0,H=0,R=0,D=0,
               ICUs=0,ICUd=0,H2=0)
    ## need to run sim for a little while to populate state!
    resICU <- run_sim(params,state,start_date=sdate,end_date="1-Apr-2020")
    state_end <- unlist(round(tail(resICU[,-1],1)))
    ## sanity check
    plot(resICU)

    ## FIXME: just populate all states
    state_end <- c(S = 19809133, E = 77785, Ia = 18177, Ip = 4699, Im = 32589, 
                   Is = 1349, H = 768, R = 55410, D = 103, ICUs = 154, ICUd = 109, 
                   H2 = 23)


    M <- make_ratemat(state_end,params2D,do_ICU=TRUE)
    M2 <- M
    M2[M2>0] <- 1 ## make all edges == 1 (otherwise edges are so narrow they disappear; could scale up or plot as flows instead if we wanted to show strength?

    image(Matrix::Matrix(M2))
    g <- igraph::graph_from_adjacency_matrix(M2)

    ## FIXME: layout??
    png("../pix/ICU_flow1.png")
    plot(g, layout=layout_as_tree)
    dev.off()
}

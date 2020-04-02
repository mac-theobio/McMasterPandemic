
##' return growth rate (from Jacobian)
##' @param p parameters
##' @export
get_r <- function(p,s) {
    res <- max(eigen(make_jac(state=s,params=p))$values)
    return(res)
}

##' return mean gen interval
##' @param p parameters
##' @export
get_GI <- function(p) {
    res <- with(as.list(p),
        1/gamma+alpha/lambda_a+
        (1-alpha)*(1/lambda_p + mu/lambda_m + (1-mu)/lambda_s))
    return(res)
}


##' compute moments of GI
##' @param params parameters
##' @export
get_giMoments <- function(params) {
    ## FIXME: assumes ICU1 model. Consider adding a test in case this changes?
    ##  (will have to rethink this once we have a structured model)
    Rv <- get_R0(params, components=TRUE)
    R <- sum(Rv)
    irates <- with(as.list(params), c(lambda_a, lambda_p, lambda_m, lambda_s))
    gamma  <- params[["gamma"]]

	 Gibar <- sum((Rv/R)/irates)
	 Givar <- sum((Rv/R)/irates^2) - Gibar^2

    Gbar <- Gibar + 1/gamma
    Gvar <- Givar + 1/gamma^2
    c(Gbar=Gbar, kappa=Gvar/Gbar^2)
}


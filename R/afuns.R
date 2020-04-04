
##' return growth rate (from Jacobian)
##' @param p parameters
##' @export
get_r <- function(p) {
    res <- max(eigen(make_jac(params=p))$values)
    return(res)
}

## OBSOLETE: naive/unweighted GI mean
## get_GI <- function(p) {
##     res <- with(as.list(p),
##         1/gamma+alpha/lambda_a+
##         (1-alpha)*(1/lambda_p + mu/lambda_m + (1-mu)/lambda_s))
##     return(res)
## }


##' compute moments of generation interval (mean and CV^2)
##' @param params parameters
##' @export
get_GI_moments <- function(params) {
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

##' calculate R0 for a given set of parameters
##' @param params parameters
##' @param components report R0 component-by-component?
##' @export
get_R0 <- function(params, components=FALSE) {
    ## FIXME: assumes ICU1 model. Consider adding a test in case this changes?
    ##  (will have to rethink this once we have a structured model)
    with(as.list(params), {
        comp <- beta0*c(alpha*Ca/lambda_a,
        (1-alpha)*c(Cp/lambda_p,mu*(1-iso_m)*Cm/lambda_m,(1-mu)*(1-iso_s)*Cs/lambda_s ))
        if (components) return(comp)
        return(sum(comp))
    })
}


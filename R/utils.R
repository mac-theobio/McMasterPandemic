
##' @importFrom methods is
## attempt convert x to a date unless it already is one
ldmy <- function(x) if (is(x,"Date")) x else lubridate::dmy(x)

## self-naming list (copied from lme4:::namedList)
nlist <- function (...) {
    L <- list(...)
    snm <- vapply(substitute(list(...)), deparse, character(1))[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

thin <- function(x,ndt=1) {
    if(ndt==1) return(x)
    x  <- x[seq(nrow(x)) %% ndt == 1,]
    return(x)
}

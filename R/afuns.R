## approximation functions

##' return growth rate (from Jacobian)
##' @param p parameters
##' @param method computation method
##' @export
get_r <- function(p, method=c("expsim","kernel","analytical")) {
    ## expsim and kernel match well, analytical is ???
    method <- match.arg(method)
    res <- switch(method,
                  analytical = {
                      warning("Jacobian-based r may be unreliable!")
                      max(eigen(make_jac(params=p))$values)
                  },
                  kernel = {
                      get_kernel_moments(p)[["r0"]]
                  },
                  expsim = {
                      rExp(p)
                  })
    return(res)
}

##' get dominant eigenvector
##' @param p parameters
##' @param method computational method
##' @export
get_evec <- function(p, method=c("expsim","analytical"),...) {
    method <- match.arg(method)
    res <- switch(method,
                  expsim=rExp(p,return_val="eigenvector",...),
                  analytical= {
                      J <- make_jac(params=p)    
                      ee <- eigen(J)
                      v <- ee$vectors
                      rownames(v) <- rownames(J)
                      dom_vec <- v[,which.max(ee$values)]
                      drop_vars <- c("date","t","S","R","D","foi","hosp","X")
                      dd <- abs(dom_vec[!names(dom_vec) %in% drop_vars])
                      dd/sum(dd)
                  })
    return(res)
}

##' compute mean generation interval from parameters
##' @param p a set of parameters
##' @param method computational method
##' @export
get_Gbar <- function(p, method=c("analytical","kernel")) {
    method <- match.arg(method)
    res <- switch(method,
                  analytical = {
                      get_GI_moments(p)[["Gbar"]]
                  },
                  kernel = {
                      ## ???
                      get_kernel_moments(p)[["Gbar"]]
                  })
    return(res)
}

##' compute R0, r, etc. based on kernel computation
##' @param params parameter vector
##' @export
get_kernel_moments <- function(params) {
    gg <- get_GI_moments(params)
    nt <- gg[["Gbar"]]*10
    kk <- transKernel(params, do_hazard=FALSE, steps=nt)$foi
    ## FIXME: check agreement between get_GI_moments() and kk ?
    return(kernelMoments(kk))
}


##' compute moments of generation interval (mean and CV^2)
##' @param params parameters
##' @export
get_GI_moments <- function(params) {
	## FIXME: assumes ICU1 model. Consider adding a test in case this changes?
	##  (will have to rethink this once we have a structured model)
	Rv <- c(0, get_R0(params, components=TRUE))
        R <- sum(Rv)
        ## FIXME: get rates, use names rather than numeric indices below
	boxtimes <- with(as.list(params), 
	 	1/c(sigma, gamma_a, gamma_p, gamma_m, gamma_s)
	)
	boxvars <- boxtimes^2
	classtimes <- c(boxtimes[[1]]
		, boxtimes[[1]] + boxtimes[[2]]
		, boxtimes[[1]] + boxtimes[[3]]
		, boxtimes[[1]] + boxtimes[[3]] + boxtimes[[4]]
		, boxtimes[[1]] + boxtimes[[3]] + boxtimes[[5]]
	)

	## Not quite right for GI!
	classvars <- c(boxvars[[1]]
		, boxvars[[1]] + boxvars[[2]]
		, boxvars[[1]] + boxvars[[3]]
		, boxvars[[1]] + boxvars[[3]] + boxvars[[4]]
		, boxvars[[1]] + boxvars[[3]] + boxvars[[5]]
	)

	Gbar <- sum((Rv/R)*classtimes)
	Gvar <- sum((Rv/R)*classvars)
return(c(R0=R, Gbar=Gbar, Gvar_approx=Gvar))
}

##' calculate R0 for a given set of parameters
##' @param params parameters
##' @param components report R0 component-by-component?
##' @param method computation method
##' @export
get_R0 <- function(params, components=FALSE,
                   method=c("analytical","kernel")) {
    method <- match.arg(method)
    res <- switch(method,
                  
	## FIXME: assumes ICU1 model. Consider adding a test in case this changes?
	##  (will have to rethink this once we have a structured model)
                  analytical = with(as.list(params), {
                      comp <- beta0*c(alpha*Ca/gamma_a,
                      (1-alpha)*c(Cp/gamma_p,
                                  mu*(1-iso_m)*Cm/gamma_m,
                                  (1-mu)*(1-iso_s)*Cs/gamma_s ))
                      if (components) comp else sum(comp)
                  }),
                  kernel=get_kernel_moments(params)[["R0"]]
                  )
    ## FIXME: helpful to name this, but we'll need something smarter in general:
    if (length(res) == 4) names(res) <- c("asymptomatic", "pre-symptomatic", "mild", "severe")
    return(res)
}

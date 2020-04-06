## kernel functions
## not exported or documented yet!

## FIXME: how to aggregate/thin/adjust time scale?
## compute kernel by brute force
transKernel <- function(par, steps=100, do_hazard=TRUE,
                        ndt=1){
        if (ndt>1) warning("ndt not fully implemented")
        par[["N"]] <- 1   ## ? redundant ?
	state <- make_state(N=1, E=1)
	return(run_sim_range(par, state
		, nt=steps*ndt
		, step_args = list(do_hazard=do_hazard)
	))
}
## allow caching of results
transKernel <- memoise::memoise(transKernel)

## FIXME: kernel should ideally be an object with k and lag
## (a class with methods)
## Building with r instead of λ for this reason

## Investigate r (combine with uniroot to get Euler's r)
discountGap <- function(r, k){
	lag <- 1:length(k)
	discountR <- sum(k*exp(-lag*r))
	return(discountR-1)
}

## Investigate κ (combine with uniroot to get κ_eff)
kappaGap <- function(kappa, rho, R){
	R_est <- (1+rho*kappa)^(1/kappa)
	return(R-R_est)
}

## moments of a kernel with parameters/convolution vector k
## FIXME: principled way to deal with lower-bound issues
kernelMoments <- function(k,lwr=0.01){
	lag <- seq(length(k))
	R0 <- sum(k)
	Gbar <- sum(k*lag)/R0
	Gvar <- sum(k*lag^2)/R0 - Gbar^2
	r0 <- (uniroot(discountGap, k=k, lower=lwr, upper=2))$root
	kappa_eff <- (uniroot(kappaGap, rho=r0*Gbar,
                              R=R0, lower=lwr, upper=2))$root

	return(c(R0=R0, Gbar=Gbar, r0=r0
		, kappa=Gvar/Gbar^2
		, kappa_eff=kappa_eff
	))
}

### What to multiply beta by to get a desired r
##' @importFrom stats uniroot
rmult <- function(k, r){
	uniroot(f=function(m) {discountGap(r, m*k)}
		, lower=1/10, upper=10
	)$root
}

## run a pure-exponential sim
## return either r or eigenvector
rExp <- function(par, steps=100, ndt=1,
                 do_hazard=FALSE,
                 return_val=c("r0","eigenvector"))
{
        return_val <- match.arg(return_val)
        if (ndt>1) warning("ndt not fully implemented")
        par[["N"]] <- 1   ## ? redundant ?
	state <- make_state(N=1, E=1e-5)
	r <- run_sim_range(par, state
                         , nt=steps*ndt
                         , step_args = list(do_hazard=do_hazard,
                                            do_exponential=TRUE))
        nn <- ndt*steps
        ## DRY: get_evec()
        drop_vars <- c("date","t","S","R","D","foi")
        uf <- function(x,pos) unlist(x[pos,!names(r) %in% drop_vars])
        r_last <- uf(r,nn)
        r_nextlast <- uf(r,nn-1)
        ret <- switch(return_val,
                      r0=mean(log(r_last/r_nextlast)),
                      eigenvector=unlist(r_last/sum(r_last)))
        return(ret)
}

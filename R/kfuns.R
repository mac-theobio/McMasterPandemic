## kernel functions
## not exported or documented yet!

## FIXME: how to aggregate/thin/adjust time scale?
## compute kernel by brute force
transKernel <- function(par, steps=100, do_hazard=TRUE,
                        ndt=1){
        if (ndt>1) warning("ndt not fully implemented")
        par[["N"]] <- 1   ## ? redundant ?
	state <- make_state(N=1, E0=1, use_eigvec=FALSE)
	return(run_sim_range(par, state
		, nt=steps*ndt
		, step_args = list(do_hazard=do_hazard)
	))
}
## allow caching of results
if (requireNamespace("memoise")) {
    transKernel <- memoise::memoise(transKernel)
}

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

##
##' run a pure-exponential sim;
##' uses run_sim_range with a population of 1 (proportions) and a very small starting value,
##'  run for 100 steps (by default)
##' used to calculate either r (technically r0) or eigenvector (for distributing initial exposed across states)
##' @param params parameters
##' @param steps number of steps to run
##' @param ndt sub-time steps
##' @param do_hazard run with hazard
##' @param testify testing compartments
##' @param return_val return growth rate or eigenvector?
##' @param type model type (passed to \code{\link{make_state}})
##' @examples
##' pp <- read_params("PHAC_testify.csv")
##' rExp(pp)
##' rExp(pp,return_val="eigenvector")
##' rExp(pp,return_val="eigenvector",testify=TRUE)
##' @export
rExp <- function(params, steps=100, ndt=1,
                 do_hazard=FALSE,
                 state = NULL,
                 testify=has_testing(params=params),
                 return_val=c("r0","eigenvector","sim"),
                 type="ICU1")
{
        return_val <- match.arg(return_val)
        if (ndt>1) warning("ndt not fully implemented")

        if(has_vax(params) && return_val == "r0"){
          if(is.null(state)) stop("must provide current state for accurate estimate of r0 with vaxified model")
          if(!isTRUE(all.equal(sum(state), 1))) stop("state vector must sum to 1")
        }

        ## need to set total population size to 1
        if(has_age(params)){
          ## with age, use a population distribution
          params[["N"]] <- mk_Nvec(attr(params, "age_cat"), Ntot = 1)
        } else {
          ## without age, just set total population size
          params[["N"]] <- 1
        }

        ## if vaxify, turn off flows between strata to keep everything constant in time
        if(has_vax(params)){
          ## turn off flows between strata
          params[["vax_doses_per_day"]] <- 0
          params[["vax_response_rate"]] <- 0
        }

        ## set up base state

        ## potential recursion here: have to make sure
        if(is.null(state)){
          state <- make_state(N=1, E0=1e-5, type="ICU1",
                              use_eigvec=FALSE,
                              params=params,
                              ageify=FALSE,
                              vaxify=FALSE,
                              testify=FALSE)  ## FIXME: don't assume ICU1?
          ## ageify and then vaxify state, as indicated by params
          if(has_age(params)){
            state <- expand_state_age(state,
                                      age_cat = attr(params, "age_cat"),
                                      Nvec = params[["N"]])
          }

          if(has_vax(params)){
            state <- expand_state_vax(state,
                                      vax_cat = attr(params, "vax_cat"))
          }
        }

      	## make initial ratemat
      	M <- make_ratemat(state=state, params=params)

      	## add testify, as indicated by params
      	if (testify) {
              M <- testify(M,params)
              state <- expand_stateval_testing(state,method="untested")
      	}
      	## run exponential simulation
      	r <- run_sim_range(
      	   params
         , state
         , nt=steps*ndt
         , M = M
         , step_args = list(do_hazard=do_hazard,
                            do_exponential=TRUE,
                            testwt_scale="none"))

      	## return whatever is being requested
        if (return_val=="sim") return(r)

      	if (return_val=="r0" && has_vax(params)) r <- condense_vax(r)

        nn <- ndt*steps
        ## DRY: get_evec()
        drop_vars <- c("date","t","D","foi","R","X","N","P","V")
        if (!testify) drop_vars <- c("S", drop_vars)
        drop_re <- paste0("^(",paste(drop_vars,collapse="|"),")")
        ## FIXME: safer version of this:
        ##   keep_vars_regexp <- "^[EIHh]"
        ##   unlist(x[pos, grepl(keep_vars_regexp, names(r))])
        uf <- function(x,pos) unlist(x[pos,!grepl(drop_re,names(r))])

        r_last <- uf(r,nn)
        r_nextlast <- uf(r,nn-1)
        r0_values <- log(r_last/r_nextlast)
        # print(r0_values)
        # print(var(r0_values))
        ## check if we've converged numerically
        ## (the discrepancy between r0_values for each compartment should be negligible, just due to numerical error)
        if(return_val =="r0" && var(r0_values)>1e-8) warning("the exponential simulation has not converged: please iterate for more steps.")
        ret <- switch(return_val,
                      ## log mean(x(t+1)/x(t))
                      r0=mean(r0_values),
                      ## normalized state vector at last time step
                      eigenvector=unlist(r_last/sum(r_last)))
        return(ret)
}

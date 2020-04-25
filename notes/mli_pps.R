## prop pred samp hacks

forecast_ensemble_mli <- function(fit,
                              nsim=200,
                              forecast_args=fit$forecast_args,
                              qvec=c(0.05,0.5,0.95),
                              qnames=c("lwr","value","upr"),
                              seed=NULL,
										equalwts = TRUE,
                              .progress=if (interactive()) "text" else "none"
                              ) {

    var <- NULL
    ## FIXME: include baseline forecast as well?
    
    ## parameter ensemble
    if (!is.null(seed)) set.seed(seed)
        ## HACK: for aggregated data, NB fit is unhappy because there is severe underdispersion (because of
        ##  overfitting to time series); NB disp parameter is >>> 1
        ##  we seem to be able to get away with ignoring it completely here
        ##  (not needed for forecast ...)
    ## might help to fix it ...
    ff <- function(p, return_val="vals_only") {
        do.call(forecast_sim,
                c(list(p,return_val=return_val),forecast_args))
    }

    ## baseline fit
    r <- ff(coef(fit$mle2), return_val="aggsim")

    ## Wald sample
    ## FIXME: count number of distribution params?
    ## FIXME: use pop_pred_samp()?
    e_pars <- as.data.frame(MASS::mvrnorm(nsim,
                                          mu=coef(fit$mle2),
                                          Sigma=bbmle::vcov(fit$mle2)))
    e_pars <- pop_pred_samp_mli(fit$mle2,n=nsim,PDify=TRUE,return_wts=TRUE)	
	 wts <- e_pars[,"wts"]
	 if(equalwts){wts <- rep(1,nsim)}
	 e_pars <- e_pars[,1:(ncol(e_pars)-1)]

    ## run for all param vals in ensemble
    ## tried with purrr::pmap but too much of a headache
    t1 <- system.time(e_res <- plyr::alply(as.matrix(e_pars)
                                         , .margins=1
                                         , .fun=ff
                                         , .progress=.progress  ))
    ## get quantiles per observation
    e_res_temp <- (e_res %>% dplyr::bind_cols())
    res_rows <- nrow(e_res_temp)
    e_res2 <- data.frame(lwr = rep(NA,res_rows)
    		, value = rep(NA,res_rows)
    		, upr = rep(NA,res_rows)
    )
    for(i in 1:res_rows){
    	if(!all(is.na(e_res_temp[i,]))){
    		e_res2[i,] <- Hmisc::wtd.quantile(as.numeric(e_res_temp[i,]),qvec,weights=wts,na.rm=TRUE)
    	}
    }
    
    ## date/var values
    e0 <- (dplyr::select(r,date,var)
        %>% dplyr::as_tibble()
    )
    ## combine quantiles with the original date/var columns
    e_res3 <- dplyr::bind_cols(e0, e_res2)
    return(e_res3)
}

pop_pred_samp_mli <- function (object, n = 1000, n_imp = n * 10, return_wts = FALSE, 
    impsamp = FALSE, PDify = FALSE, PDmethod = NULL, tol = 1e-06, 
    return_all = FALSE, rmvnorm_method = c("mvtnorm", "MASS"), 
    fix_params = NULL) 
{
    rmvnorm_method <- match.arg(rmvnorm_method)
    min_eval <- function(x) {
        ev <- eigen(x, only.values = TRUE)$values
        if (is.complex(ev)) {
            print(x)
            print(ev)
            stop("covariance matrix with complex eigenvalues (!)")
        }
        min(ev)
    }
    cc_full <- object@fullcoef
    cc <- object@coef
    keep_params <- !names(cc) %in% fix_params
    cc <- cc[keep_params]
    vv <- vcov(object)
    vv <- vv[keep_params, keep_params]
    Lfun <- object@minuslogl
    fixed_pars <- setdiff(names(object@fullcoef), names(cc))
    res <- matrix(NA, nrow = n, ncol = length(cc_full), dimnames = list(NULL, 
        names(cc_full)))
    if (any(is.na(cc))) 
        return(res)
    bad_vcov <- any(is.na(vv))
    if (!bad_vcov) {
        min_eig <- min_eval(vv)
    }
    else {
        min_eig <- NA
    }
    if (is.na(min_eig) || any(min_eig < tol)) {
        if (!PDify) {
            stop("NA or non-positive definitive variance-covariance matrix ", 
                sprintf("(min eig=%f): ", min_eig), "consider PDify=TRUE (and probably impsamp=TRUE)")
        }
        hh <- object@details$hessian[keep_params, keep_params]
        if (any(is.na(hh))) {
            warning("NA values in Hessian set to zero: check results *very* carefully!")
            hh[is.na(hh)] <- 0
        }
        if ((is.null(PDmethod) && bad_vcov) || identical(PDmethod, 
            "King")) {
            warning("using EXPERIMENTAL King et al method")
            vv <- crossprod(as.matrix(bdsmatrix::gchol(MASS::ginv(hh)), 
                ones = FALSE))
        }
        else {
            vv <- as.matrix(Matrix::nearPD(vv)$mat)
        }
    }
    mv_n <- if (impsamp) 
        n_imp
    else n
    res[, names(cc)] <- mv_vals <- switch(rmvnorm_method, mvtnorm = mvtnorm::rmvnorm(mv_n, 
        mean = cc, sigma = vv), MASS = MASS::mvrnorm(mv_n, mu = cc, 
        Sigma = vv))
    if (length(fixed_pars) > 0) {
        for (p in fixed_pars) {
            res[, p] <- object@fullcoef[p]
        }
    }
    if (!(impsamp || return_wts)) 
        return(res)
    mv_wts <- mvtnorm::dmvnorm(mv_vals, mean = cc, sigma = vv, 
        log = TRUE)
    if (all(is.na(mv_wts)) && length(mv_wts) == 1) {
        mv_wts <- rep(NA, length(mv_vals))
        warning("can't compute MV sampling probabilities")
    }
    # L_wts0 <- -1 * apply(res, 1, Lfun)
    L_wts0 <- -1 * apply(res,1,function(x)Lfun(p=x
		, opt_pars=opt_pars
		, base_params=params
		, start_date=start_date
		, end_date = end_date
		, data=object@data$data)
    )

    L_wts <- L_wts0 - mv_wts
    L_wts <- exp(L_wts - max(L_wts, na.rm = TRUE))
    L_wts <- L_wts/sum(L_wts, na.rm = TRUE)
    eff_samp <- 1/sum(L_wts^2, na.rm = TRUE)
    res <- cbind(res, wts = L_wts)
    attr(res, "eff_samp") <- eff_samp
    if (return_all) {
        return(cbind(res, loglik = L_wts0, mvnloglik = mv_wts))
    }
    if (return_wts) 
        return(res)
    res <- res[sample(seq(nrow(res)), size = n, prob = L_wts, 
        replace = TRUE), ]
    return(res)
}

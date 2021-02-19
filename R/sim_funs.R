##' construct Jacobian matrix for ICU model
##' (not quite complete: doesn't include flows to R)
## FIXME: derive from make_ratemat 
##' @param state state vector (named)
##' @param params parameter vector
##' @export
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params[["N"]],E0=params[["E0"]])
##' ## state[c("E","Ia","Ip","Im","Is")] <- 1
##' state[["E"]] <- 1
##' J <- make_jac(params,state)
##' J["S","S"]
##' Jr <- J[1:6,1:6]
##' round(Jr,3)
##' eigen(Jr)$values
##' make_jac(params)
make_jac <- function(params, state=NULL) {
    ## circumvent test code analyzers ... problematic ...
    S <- E <- Ia <- Ip <- Im <- Is <- H <- NULL
    H2 <- ICUs <- ICUd <- D <- R <- beta0 <- Ca <- Cp  <- NULL
    Cm <- Cs <- alpha <- sigma <- gamma_a <- gamma_m <- gamma_s <- gamma_p  <- NULL
    rho <- delta <- mu <- N <- E0 <- iso_m <- iso_s <- phi1  <- NULL
    phi2 <- psi1 <- psi2 <- psi3 <- c_prop <- c_delaymean <- c_delayCV  <- NULL
    ####
    if (is.null(state)) {
        state <- make_state(N=params[["N"]],E0=1e-3,
                            use_eigvec=FALSE)
    }
    np <- length(params)
    ns <- length(state)
    ## make state and param names locally available (similar to with())
    P <- c(as.list(state),as.list(params))
    unpack(P) ## extract variables
    ## blank matrix
    M <- matrix(0,
                nrow=ns, ncol=ns,
                dimnames=list(from=names(state),to=names(state)))
    Ivec <- c(Ia, Ip, Im,Is)
    Iwt <- beta0/N*c(Ia=Ca,Ip=Cp,Im=(1-iso_m)*Cm,Is=(1-iso_s)*Cs)
    Ivars <- c("Ia","Ip","Im","Is")
    M["S","S"] <- -sum(Ivec*Iwt)
    M["S",Ivars] <- -S*Iwt[Ivars]
    M["E",c("S",Ivars)] <- -M["S",c("S",Ivars)]
    M["E","E"] <- -sigma
    M["Ia","E"] <- alpha*sigma
    M["Ia","Ia"] <- -gamma_a
    M["Ip","E"] <- (1-alpha)*sigma
    M["Ip","Ip"] <- -gamma_p
    M["Im","Ip"] <- mu*gamma_p
    M["Im","Im"] <- -gamma_m
    M["Is","Ip"] <- (1-mu)*gamma_p
    M["Is","Is"] <- -gamma_s
    M["H","Is"]  <- phi1*gamma_s
    M["H","H"]   <- -rho
    M["ICUs","Is"] <- (1-phi1)*(1-phi2)*gamma_s
    M["ICUs","ICUs"] <- -psi1
    M["H2","ICUs"] <- psi1
    M["H2","H2"] <- -psi3
    M["ICUd","Is"] <- (1-phi1)*phi2*gamma_s
    M["ICUd","ICUd"] <- -psi2
    M["D","ICUd"] <- psi2
    M["R","Ia"] <- gamma_a
    M["R","Im"] <- gamma_m
    M["R","H"] <- rho
    M["R","H2"] <- psi3
    return(M)
}

##' construct vector of transmission multipliers
##' @param state state vector
##' @param params parameter vector
##' @param full include non-infectious compartments (with transmission of 0) as well as infectious compartments?
##' @export
## QUESTION: is the main testify argument to this function used?
make_betavec <- function(state, params, full=TRUE) {
    Icats <- c("Ia","Ip","Im","Is")
    testcats <- c("_u","_p","_n","_t")
    ## NB meaning of iso_* has switched from Stanford model
    ## beta_vec0 is the vector of transmission parameters that apply to infectious categories only
    beta_vec0 <- with(as.list(params),
                      beta0/N*c(Ca,Cp,(1-iso_m)*Cm,(1-iso_s)*Cs))
    names(beta_vec0) <- Icats
    if (has_age(params)) {
        Cmat <- params$Cmat
        a_names <- rownames(Cmat)
        new_names <- expand_names(Icats, a_names)
        beta_vec0 <- t(kronecker(Cmat,matrix(beta_vec0)))
        dimnames(beta_vec0) <- list(a_names,new_names)
        beta_vec0 <- Matrix(beta_vec0)
    }
    ## assume that any matching values will be of the form "^%s_" where %s is something in Icats
    ## lapply(Icats, function(x) grep(sprintf("^%s_"), names(state))
    ## FIXME: we should be doing this by name, not assuming that all infectious compartments are expanded
    ##  into exactly 4 subcompartments, in order (but this should work for now??)
    if (has_testing(state=state)) {  ## testified!
        if (has_age(params)) stop("can't combine age and testing yet")
        beta_vec0 <- rep(beta_vec0,each=length(testcats))
        names(beta_vec0) <- unlist(lapply(Icats,function(x) paste0(x,testcats)))
        ## FIXME: also adjust _n, _p components?
        pos_vals <- grep("_t$",names(beta_vec0))
        beta_vec0[pos_vals] <- beta_vec0[pos_vals]*(1-params[["iso_t"]])
    }
    if (!full) return(beta_vec0)
    ## By default, make a vector of zeroes for all the states,
    ## then fill in infectious ones
    if (!has_age(params)) {
        beta_vec <- setNames(numeric(length(state)),names(state))
        beta_vec[names(beta_vec0)] <- beta_vec0
    } else {
        beta_vec <- matrix(0, nrow=nrow(beta_vec0), ncol=length(state),
                           dimnames=list(rownames(beta_vec0), names(state)))
        beta_vec[rownames(beta_vec0),colnames(beta_vec0)] <- matrix(beta_vec0)
    }
    return(beta_vec)
}

## make_ratemat()
##' Create transition matrix
##'
##' Defines rates (per day) of flow \emph{from} compartment \code{i}
##' (row) \emph{to} compartment \code{j} (column).
##'
##' @details
##' The rates are as follows:
##' 
##' \eqn{ S to E:  - (\beta_0 / N) S (C_a I_a + C_p I_p + (1-iso_m)C_m I_m + (1-iso_s)C_s I_s) }
##' \eqn{ E to I_a: }
##' \eqn{ E to I_p: }
##' \eqn{ ... }
##'
##' See \code{\link{read_params}} for parameter definitions.
##' 
##' @note
##' Base version matches structure of Stanford/Georgia models
##' \itemize{
##'   \item flow diagram: see \url{http://covid-measures.stanford.edu/} 'model details' tab
##'         or \code{../pix/model_schematic.png}
##'   \item parameter definitions: see \code{params_CI_base.csv}, \code{params_ICU_diffs.csv}
##' }
##' 
##' @param state named vector of states
##' @param params named vector of parameters
##' @param do_ICU include additional health utilization compartments
##' @param sparse return sparse matrix?
##' @param symbols return character (symbol) form? (FIXME: call adjust_symbols here rather than in show_ratemat()?)
##' @return matrix of (daily) transition rates
##  *need* Matrix version of rowSums imported to handle sparse stuff below!! 
##' @importFrom Matrix Matrix rowSums colSums
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params[["N"]],E0=params[["E0"]])
##' M <- make_ratemat(state,params)
##' if (require(Matrix)) {
##'    image(Matrix(M))
##' }
##' make_ratemat(state,params,symbols=TRUE)
##' @export
make_ratemat <- function(state, params, do_ICU=TRUE, sparse=FALSE,
                         symbols=FALSE) {
    ## circumvent test code analyzers ... problematic ...
    alpha <- sigma <- gamma_a <- gamma_m <- gamma_s <- gamma_p  <- NULL
    rho <- delta <- mu <- N <- E0 <- iso_m <- iso_s <- phi1  <- NULL
    phi2 <- psi1 <- psi2 <- psi3 <- c_prop <- c_delaymean <- c_delayCV  <- NULL
    ## default values, will be masked (on purpose) by unpacking params/state
    nonhosp_mort <- 0
    ####
    np <- length(params)
    if (is.list(params)) {
        nps <- lengths(params)
        if (has_age(params)) {
            na <- nrow(params[["Cmat"]])
            bad_len <- which(!nps %in% c(1,na,na^2))
            if (length(bad_len)>0) {
                stop(sprintf("elements of params must be length 1, %d or %d: %s",
                             na,na^2,
                             paste(names(params)[bad_len],collapse=", ")))
            }
        } else {
            ## FIXME: better error message ...
            stopifnot(all(nps==1))
        }
    }
    state_names <- untestify_statenames(names(state))
    ns <- length(state_names)
    ## make param names locally available (similar to with())
    ## DON'T unpack states, we don't need them
    ## (the only state-dependent per capita rates are testing
    ## and infection, those get handled elsewhere)
    P <- as.list(params) 
    unpack(P)
    ## blank matrix
    M <- matrix(0,
                nrow=ns, ncol=ns,
                dimnames=list(from=state_names, to=state_names))

    ## generic assignment function, indexes by regexp rather than string
    afun <- function(from, to, val) {
        if (!symbols) {
            M[pfun(from, to, M)] <<- val
        } else {
            M[pfun(from,to, M)] <<- deparse(substitute(val))            
        }
    }
    
    ## fill entries
    beta_vec <- make_betavec(state,params)
    ## FIXME: call update_foi() here?
    if (!has_age(params)) {
        afun("S", "E", sum(beta_vec*state[names(beta_vec)]))
    } else {
        afun("S", "E", beta_vec %*% state[colnames(beta_vec)])
    }
    afun("E", "Ia", alpha*sigma)
    afun("E","Ip", (1-alpha)*sigma)
    afun("Ia","R", gamma_a)
    afun("Ip","Im", mu*gamma_p)
    afun("Ip","Is", (1-mu)*gamma_p)
    afun("Im","R", gamma_m)
    if (!do_ICU) {
        ## simple hospital model as in Stanford/CEID
        afun("Is","H", gamma_s)
        afun("H","D", delta*rho)
        afun("H","R", (1-delta)*rho)
    } else {
        ## FIXME: A better term than "acute" to mean the opposite of intensive?
        ## four-way split (direct to D, acute care, ICU/survive, ICUD/die)?
        afun("Is","H", (1-nonhosp_mort)*phi1*gamma_s)
        afun("Is","ICUs", (1-nonhosp_mort)*(1-phi1)*(1-phi2)*gamma_s)
        afun("Is","ICUd", (1-nonhosp_mort)*(1-phi1)*phi2*gamma_s)
        afun("Is","D", nonhosp_mort*gamma_s)
        afun("ICUs","H2", psi1) ## ICU to post-ICU acute care
        afun("ICUd","D", psi2)  ## ICU to death
        afun("H2","R", psi3)  ## post-ICU to discharge
        ## H now means 'acute care' only; all H survive & are discharged
        afun("H","D", 0)
        afun("H","R", rho) ## all acute-care survive
        if (any(grepl("^X",colnames(M)))) {
            ## FIXME: check for age?
            afun("Is","X", M[pfun("Is","H",M)]) ## assuming that hosp admissions mean *all* (acute-care + ICU)
        }
    }
    if (sparse) {
        M <- Matrix::Matrix(M)
    }
    return(M)
}

##' calculate only updated force of infection
##' at present, this is the only state-dependent \emph{per capita} rate
##' maybe more efficient than modifying & returning the whole matrix
##' @inheritParams make_ratemat
##' @param beta_vec vector of transmission rates (matching state vector)
##' @export
## FIXME DRY from make_ratemat
update_foi <- function(state, params, beta_vec) {
    ## update infection rate
    if (is.matrix(beta_vec)) {
        ## FIXME, check dimensions etc.
        foi <- beta_vec %*% state
    } else {
        if(length(state) != length(beta_vec)){
            stop("length of state and beta_vec are not the same")
        }
        foi <- sum(state*beta_vec)
    }
    if (has_zeta(params)) {
        ## suppose zeta is a vector zeta1, zeta2, zeta3, ...
        ##  we also need parameters   zeta_break1, zeta_break2 (0<zbi<1)
        ##  one *fewer* break parameter than zeta_i value
        ## if 0< S/N < zeta_break1   -> zeta1
        ##  zeta_break1 < S/N < zeta_break2 -> zeta2
        ## ...
        ##  zeta_breakx < S/N < 1  -> zetax
        Susc_frac <- 1/params[["N"]]*sum(state[grep("^S_?",names(state))])
        if (any(grepl("zeta[0-9]",names(params)))) {
            zeta <- with(as.list(params),
                         if (Susc_frac<zeta_break) zeta1 else zeta2)
        }
        ## alternately could just make it a vector
        ## ... but this messes with age-structured stuff
        ## if (length(zeta)>0) {
        ## ...
        ## }
        ## ???? het at pop level or sub-category level??
        foi <- foi*with(as.list(params), Susc_frac^zeta)
    }
    return(foi)
}

update_ratemat <- function(ratemat, state, params, testwt_scale="N") {
    if (inherits(ratemat,"Matrix")) {
        aa <- c("wtsvec","posvec","testing_time")
        saved_attrs <- setNames(lapply(aa,attr,x=ratemat),aa)
    }
    ## update testing flows. DO THIS FIRST, before updating foi: **** assignment via cbind() to Matrix objects loses attributes???
    if (has_testing(state)) {
        testing_time <- attr(ratemat,"testing_time")  ## ugh  (see **** above)
        ## positions of untested, positive-waiting, negative-waiting compartments
        ## (flows from _n, _p to _t, or back to _u, are always at per capita rate omega, don't need
        ##  to be updated for changes in state)
        ## FIXME: backport to testify?
        u_pos <- grep("_u$",rownames(ratemat))
        p_pos <- grep("_p$",rownames(ratemat))
        n_pos <- grep("_n$",rownames(ratemat))
        ## original/unscaled prob of positive test by compartment, testing weights by compartment
        posvec <- attr(ratemat,"posvec")
        if (is.null(posvec)) stop("expected ratemat to have a posvec attribute")
        wtsvec <- attr(ratemat,"wtsvec")
        if (is.null(wtsvec)) stop("expected ratemat to have a wtsvec attribute")
        ## scaling ...
        testing_intensity <- params[["testing_intensity"]]
        testing_tau <- params[["testing_tau"]]
        N0 <- params[["N"]]
        W <- sum(wtsvec*state[u_pos])
        sc <- switch(testwt_scale,
                     none=1,
                     N=N0/W,
                     sum_u=sum(state[u_pos])/W,
                     sum_smooth={
                         rho <- testing_intensity
                         tau <- testing_tau
                         tau*N0/(tau*W + rho*N0)
			 ## NOTE 'smoothing' doc has numerator rho*tau*N0,
                         ## but testing intensity (rho) is included in ratemat 
                         ## calculation below ...
                     })
        ratemat[cbind(u_pos,n_pos)] <- testing_intensity*sc*wtsvec*(1-posvec)
        ratemat[cbind(u_pos,p_pos)] <- testing_intensity*sc*wtsvec*posvec
        if (testing_time=="sample") {
            N_pos <- which(rownames(ratemat)=="N")
            P_pos <- which(rownames(ratemat)=="P")
            ratemat[cbind(u_pos,N_pos)] <- ratemat[cbind(u_pos,n_pos)]
            ratemat[cbind(u_pos,P_pos)] <- ratemat[cbind(u_pos,p_pos)]
        }
    }
    ratemat[pfun("S","E",ratemat)]  <- update_foi(state,params,make_betavec(state,params))
    ## ugh, restore attributes if necessary
    if (inherits(ratemat,"Matrix")) {
        for (a in aa) {
            attr(ratemat,a) <- saved_attrs[[a]]
        }
    }
    return(ratemat)
}

## do_step()
##' Take a single simulation time step
##' @inheritParams make_ratemat
##' @param ratemat transition matrix
##' @param dt time step (days)
##' @param do_hazard use hazard calculation?
##' @param do_exponential prevent outflow of susceptibles, to create a pure-exponential process?
##' @param stoch_proc stochastic process error?
##' @param testwt_scale how to scale testing weights? "none"=use original weights as specified;
##' "N" = multiply by (pop size)/(sum(wts*state[u_pop])); "sum_u" = multiply by (sum(state[u_pop])/(sum(wts*state[u_pop])))
##' @export
##' @examples
##' params1 <- read_params("ICU1.csv")
##' state1 <- make_state(params=params1)
##' M <- make_ratemat(params=params1, state=state1)
##' s1A <- do_step(state1,params1, M, stoch_proc=TRUE)
do_step <- function(state, params, ratemat, dt=1,
                    do_hazard=FALSE, stoch_proc=FALSE,
                    do_exponential=FALSE,
                    testwt_scale="N") {
    
    x_states <- c("X","N","P")                  ## weird parallel accumulators
    p_states <- exclude_states(names(state), x_states)
    ## FIXME: check (here or elsewhere) for non-integer state and process stoch?
    ## cat("do_step beta0",params[["beta0"]],"\n")
    ratemat <- update_ratemat(ratemat, state, params, testwt_scale=testwt_scale)
    if (!stoch_proc || (!is.null(s <- params[["proc_disp"]]) && s<0)) {
        if (!do_hazard) {
            ## from per capita rates to absolute changes
            flows <- sweep(ratemat, state, MARGIN=1, FUN="*")*dt
        } else {
            ## FIXME: change var names? {S,E} is a little confusing (sum, exp not susc/exposed)
            ## use hazard function: assumes exponential change
            ## (constant per capita flows) rather than linear change
            ## (constant absolute flows) within time steps
            ## handle this as in pomp::reulermultinom,
            ## i.e.
            ##    S = sum(r_i)   ## total rate
            ##    p_{ij}=(1-exp(-S*dt))*r_j/S
            ##    p_{ii}= exp(-S*dt)
            S <- rowSums(ratemat)
            E <- exp(-S*dt)
            ## prevent division-by-0 (boxes with no outflow) problems (FIXME: DOUBLE-CHECK)
            norm_sum <- ifelse(S==0, 0, state/S)
            flows <- (1-E)*sweep(ratemat, norm_sum, MARGIN=1, FUN="*")
            diag(flows) <- 0  ## no flow
        }
    } else {
        flows <- ratemat ## copy structure
        flows[] <- 0
        for (i in seq(length(state))) {
            ## FIXME: allow Dirichlet-multinomial ?
            dW <- dt
            if (!is.na(proc_disp <- params[["proc_disp"]])) {
                dW <- pomp::rgammawn(sigma=proc_disp,dt=dt)
            }
            ## FIXME: need to adjust for non-conserving accumulators!
            flows[i,-i] <- pomp::reulermultinom(n=1,
                                 size=state[[i]],
                                 rate=ratemat[i,-i],
                                 dt=dW)
            
        }
    }

    if (!do_exponential) {
        outflow <- rowSums(flows[,p_states])
    } else {
        ## want to zero out outflows from S to non-S compartments
        ##  (but leave the inflows - thus we can't just zero these flows
        ##   out in the rate matrix!)
        S_pos <- grep("^S",rownames(ratemat), value=TRUE)
        notS_pos <- grep("^[^S]",colnames(ratemat), value=TRUE)
        notS_pos <- setdiff(notS_pos, x_states)
        outflow <- setNames(numeric(ncol(flows)),colnames(flows))
        ## only count outflows from S_pos to other S_pos (e.g. testing flows)
        outflow[S_pos] <- rowSums(flows[S_pos,S_pos,drop=FALSE])
        ## count flows from infected etc. to p_states (i.e. states that are *not* parallel accumulators)
        outflow[notS_pos] <- rowSums(flows[notS_pos,p_states])
    }
    inflow <-  colSums(flows)
    state <- state - outflow + inflow
    ## check conservation (*don't* check if we are doing an exponential sim, where we
    ##  allow infecteds to increase without depleting S ...)
    MP_badsum_action <- getOption("MP_badsum_action","warning")
    MP_badsum_tol <- getOption("MP_badsum_tol",1e-12)
    if (!do_exponential
        && !(MP_badsum_action=="ignore")
        && !stoch_proc)    ## temporary: adjust reulermultinom to allow for x_states ...
    {
        calc_N <- sum(state[p_states])
        if (!isTRUE(all.equal(calc_N,params[["N"]], tolerance=MP_badsum_tol))) {
            msg <- sprintf("sum(states) != original N (delta=%1.2g)",params[["N"]]-calc_N)
            get(MP_badsum_action)(msg)
        }
    } ## not exponential run or stoch proc or ignore-sum
    return(state)
}

## run_sim()
##' Run pandemic simulation
##' @inheritParams do_step
##' @param start_date starting date (Date or character, any sensible D-M-Y format)
##' @param end_date ending date (ditto)
##' @param params_timevar three-column data frame containing columns 'Date'; 'Symbol' (parameter name/symbol); 'Relative_value' (value \emph{relative to baseline})
##' @param dt time step for \code{\link{do_step}}
##' @param ratemat_args additional arguments to pass to \code{\link{make_ratemat}}
##' @param step_args additional arguments to pass to \code{\link{do_step}}
##' @param ndt number of internal time steps per time step
##' @param stoch a logical vector with elements "obs" (add obs error?) and "proc" (add process noise?)
##' @param stoch_start dates on which to enable stochasticity (vector of dates with names 'proc' and 'obs')
##' @param condense if \code{TRUE}, use \code{\link{condense.pansim}} to reduce the number of variables in the output (in particular, collapse subclasses and return only one \code{I}, \code{H}, and \code{ICU} variable)
##' @param condense_args arguments to pass to \code{\link{condense}} (before adding observation error)
##' @param use_ode integrate via ODE rather than discrete step?
##' @param ode_args additional arguments to \code{\link[deSolve]{ode}}
##' @examples
##' params <- read_params("ICU1.csv")
##' paramsS <- update(params,c(proc_disp=0.1,obs_disp=100))
##' paramsSz <- update(paramsS, zeta=5)
##' state <- make_state(params=params)
##' time_pars <- data.frame(Date=c("2020-Mar-20","2020-Mar-25"),
##'                        Symbol=c("beta0","beta0"),
##'                        Relative_value=c(0.7,0.1),
##'                        stringsAsFactors=FALSE)
##' res1 <- run_sim(params,state,start_date="2020-Feb-1",end_date="2020-Jun-1")
##' res1X <- run_sim(params,state,start_date="2020-Feb-1",end_date="2020-Jun-1",
##'                  condense_args=list(keep_all=TRUE))
##' res1_S <- update(res1, params=paramsS, stoch=c(obs=TRUE, proc=TRUE))
##' res1_t <- update(res1, params_timevar=time_pars)
##' res1_S_t <- update(res1_S, params_timevar=time_pars)
##' res2_S_t <- update(res1_S_t,params=update(paramsS, proc_disp=0.5))
##' res3_S_t <- update(res2_S_t,stoch_start="2020-Apr-1")
##' res3_Sz <- update(res1_S, params=paramsSz)
##' plot(res3_Sz,log=TRUE,log_lwr=1e-4)
##' @importFrom stats rnbinom na.exclude napredict
##' @importFrom anytime anydate
##' @param verbose print messages (e.g. about time-varying parameters)?
##' @export
## FIXME: params_timevar
##   change param name to something less clunky? case-insensitive/partial-match columns? allow Value and Relative_value? (translate to one or the other at R code level, for future low-level code?)
## FIXME: automate state construction better
run_sim <- function(params
        , state=NULL
        , start_date="2020-Mar-20"
        , end_date="2020-May-1"
        , params_timevar=NULL
        , dt=1
        , ndt=1  ## FIXME: change default after testing?
        , stoch=c(obs=FALSE,proc=FALSE)
        , stoch_start=NULL
        , ratemat_args=NULL
        , step_args=list()
        , ode_args = list()
        , use_ode = FALSE
        , condense = TRUE
        , condense_args=NULL
        , verbose = FALSE
) {
    call <- match.call()
    
    if (is.na(params[["N"]])) stop("no population size specified; set params[['N']]")
    ## FIXME: *_args approach (specifying arguments to pass through to
    ##  make_ratemat() and do_step) avoids cluttering the argument
    ##  list, but may be harder to translate to lower-level code
    if (dt!=1) warning("nothing has been tested with dt!=1")
    start_date <- anydate(start_date); end_date <- anydate(end_date)
    if (!is.null(stoch_start)) {
        ## anydate() strips names ...
        stoch_start <- setNames(anydate(stoch_start),names(stoch_start))
    }
    if (length(stoch_start)==1) stoch_start <- c(obs=stoch_start, proc=stoch_start)
    date_vec <- seq(start_date,end_date,by=dt)
    nt <- length(date_vec)
    step_args <- c(step_args, list(stoch_proc=stoch[["proc"]]))
    drop_last <- function(x) { x[seq(nrow(x)-1),] }
    if (is.null(state)) state <- make_state(params=params, testify=FALSE)
    M <- make_ratemat(state=state, params=params)
    if (has_testing(params=params)) {
        if (!is.null(ratemat_args$testify)) {
            warning("'testify' no longer needs to be passed in ratemat_args")
        }
        testing_time <- ratemat_args$testing_time
        if (is.null(testing_time)) {
            warning("setting testing time to 'sample'")
            testing_time <- 'sample'
        }
        M <- testify(M,params,testing_time=testing_time)
        state <- expand_stateval_testing(state, params=params)
    }
    state0 <- state
    params0 <- params ## save baseline (time-0) values
    ## no explicit switches, and (no process error) or (process error for full time);
    ## we will be able to run the whole sim directly
    if (is.null(params_timevar) && (!stoch[["proc"]] || is.null(stoch_start))) {
        switch_times <- NULL
    } else {
        if (is.null(params_timevar)) {
            ## starting times for process/obs error specified, but no other time-varying parameters;
            ##  we need an empty data frame with the right structure so we can append the process-error switch times
            params_timevar <- dfs(Date=as.Date(character(0)),Symbol=character(0),Relative_value=numeric(0))
        } else {
            ## check column names
            ## FIXME:: use case-insensitive matching (apply tolower() throughout) to allow some slop?
            npt <- names(params_timevar)
            if (!all(c("Date","Symbol","Relative_value") %in% npt)) {
                stop("bad names in params_timevar: ",paste(npt,collapse=","))
            }
            params_timevar$Date <- anydate(params_timevar$Date)
            ## tryCatch(
            ##     params_timevar$Date <- as.Date(params_timevar$Date),
            ##     error=function(e) stop("Date column of params_timevar must be a Date, or convertible via as.Date"))
            params_timevar <- params_timevar[order(params_timevar$Date),]
        }
        ## append process-observation switch to timevar
        if (stoch[["proc"]] && !is.null(stoch_start)) {
            params_timevar <- rbind(params_timevar,
                                    dfs(Date=stoch_start[["proc"]],
                                        Symbol="proc_disp",Relative_value=1))
            params[["proc_disp"]] <- -1 ## special value: signal no proc error
        }
        switch_dates <- params_timevar[["Date"]]
        ## match specified times with time sequence
        switch_times <- match(switch_dates, date_vec)
        if (any(is.na(switch_times))) {
            bad <- which(is.na(switch_times))
            stop("non-matching dates in params_timevar: ",paste(switch_dates[bad], collapse=","))
        }
        if (any(switch_times==length(date_vec))) {
            ## drop switch times on final day
            warning("dropped switch times on final day")
            switch_times <- switch_times[switch_times<length(date_vec)]
        }
    } ## steps

    if (is.null(switch_times)) {
        res <- thin(ndt=ndt,
                    do.call(run_sim_range,
                            nlist(params
                                , state
                                , nt=nt*ndt
                                , dt=dt/ndt
                                , M
                                , use_ode
                                , ratemat_args
                                , step_args
                            )))
    } else {
        t_cur <- 1
        ## want to *include* end date 
        switch_times <- switch_times + 1
        ## add beginning and ending time
        times <- c(1,unique(switch_times),nt+1)
        resList <- list()
        for (i in seq(length(times)-1)) {
            for (j in which(switch_times==times[i])) {
                ## reset all changing params
                s <- params_timevar[j,"Symbol"]
                v <- params_timevar[j,"Relative_value"]
                params[[s]] <- params0[[s]]*v
                if (s=="proc_disp") {
                    state <- round(state)
                }
                if (verbose) cat(sprintf("changing value of %s from original %f to %f at time step %d\n",
                                         s,params0[[s]],params[[s]],i))
                
                ## FIXME: so far still assuming that params only change foi
                ## if we change another parameter we will have to recompute M 
            }

            resList[[i]] <- drop_last(
                thin(ndt=ndt,
                     do.call(run_sim_range,
                        nlist(params
                            , state
                            , nt=(times[i+1]-times[i]+1)*ndt
                            , dt=dt/ndt
                            , M
                            , use_ode
                            , ratemat_args
                            , step_args
                            , ode_args
                              )))
            )
            state <- attr(resList[[i]],"state")
            t_cur <- times[i]
        }
        ## combine
        res <- do.call(rbind,resList)
        ## add last row
        ## res <- rbind(res, attr(resList[[length(resList)]],"state"))
    }
    ## drop internal stuff
    ## res <- res[,setdiff(names(res),c("t","foi"))]
    res <- dfs(date=seq(start_date,end_date,by=dt),res)
    res <- res[,!names(res) %in% "t"]  ## we never want the internal time vector
    ## condense here
    if (condense) {
        res <- do.call(condense.pansim,c(list(res,params=params0,
                                              cum_reports=FALSE,
                                              het_S=has_zeta(params0)),
                                         condense_args))
    }
    if (stoch[["obs"]]) {
        if (has_zeta(params)) params[["obs_disp_hetS"]] <- NA  ## hard-code skipping obs noise
        ## do observation error here
        ## FIXME: warn on mu<0 ? (annoying with ESS machinery)
        m <- res[,-1]   ## drop time component
        if (!is.null(stoch_start)) {
            ## only add stochastic obs error after specified date
            m <- m[res$date>stoch_start[["obs"]],]
        }
        m_rows <- nrow(m)
        for (i in seq(ncol(m))) {
            nm <- names(m)[i]
            if ((vn <- paste0("obs_disp_",nm)) %in% names(params)) {
                ## variable-specific dispersion specified
                d <- params[[vn]]
            } else d <- params[["obs_disp"]]
            if (!is.na(d)) {
                ## rnbinom, skipping NA values (as below)
                m[[i]] <- suppressWarnings(rnbinom(m_rows, mu=m[[i]], size=d))
            }
        }
        res[seq(nrow(res)-m_rows+1,nrow(res)),-1] <- m
    }
    ## add cum reports *after* adding obs error
    if ("report" %in% names(res)) res$cumRep <- cumsum(ifelse(!is.na(res$report), res$report, 0))
    if ("death" %in% names(res)) res$D <- cumsum(ifelse(!is.na(res$death), res$death, 0))
    ## store everything as attributes
    attr(res,"params") <- params0
    attr(res,"state0") <- state0
    attr(res,"start_date") <- start_date
    attr(res,"end_date") <- end_date
    attr(res,"call") <- call
    attr(res,"params_timevar") <- params_timevar
	 ## attr(res,"final_state") <- state
    class(res) <- c("pansim","data.frame")
    return(res)
}



##' retrieve parameters from a CSV file
##'
##' @details
##' The parameters that must be set are:
##'
##' \eqn{   N:  }  population size
##' 
##' \eqn{   \beta_0:  }  transmission rate
##' 
##' \eqn{   1/\sigma:  }  mean \emph{latent} period
##' 
##' \eqn{   1/\gamma_a:  }  mean \emph{infectious} period for asymptomatic individuals
##' 
##' \eqn{   ... }
##' 

##' generate initial state vector
##' @param N population size
##' @param E0 initial number exposed
##' @param type (character) specify what model type this is intended
##'     for (e.g., \code{"ICU1"}, \code{"CI"}); determines state names
##' @param state_names vector of state names, must include S and E
##' @param params parameter vector (looked in for N and E0)
##' @param x proposed (named) state vector; missing values will be set
##' @param use_eigvec use dominant eigenvector to distribute non-Susc values
##'     to zero
##' @param testify expand state vector to include testing compartments (untested, neg waiting, pos waiting, pos received) ?
##' @note \code{"CI"} refers to the Stanford group's
##'     "covid intervention" model.
##' @export
##' @examples
##' p <- read_params("ICU1.csv")
##' make_state(N=1e6,E0=1)
##' make_state(params=p)
##  FIXME: can pass x, have a name check, fill in zero values
make_state <- function(N=params[["N"]],
                       E0=params[["E0"]],
                       type="ICU1h",
                       state_names=NULL,
                       use_eigvec=TRUE,
                       params=NULL,
                       x=NULL,
                       testify=NULL) {
    if (is.null(testify)) testify <- !is.null(params) && has_testing(params=params)
    if (use_eigvec && is.null(params)) stop("must specify params")
    ## select vector of state names
    state_names <- switch(type,
                          ## "X" is a hospital-accumulator compartment (diff(X) -> hosp)
                          ICU1h = c("S","E","Ia","Ip","Im","Is","H","H2","ICUs","ICUd", "D","R","X"),
                          ICU1 = c("S","E","Ia","Ip","Im","Is","H","H2","ICUs","ICUd", "D","R"),
                          CI =   c("S","E","Ia","Ip","Im","Is","H","D","R"),
                          stop("unknown type")
                          )
    state <- setNames(numeric(length(state_names)),state_names)
    if (testify) state <- expand_stateval_testing(state,method="untested")
    if (is.null(x)) {
        ## state[["S"]] <- round(N-E0)
        if (!use_eigvec) {
            ## set **first** E compartment
            state[[grep("E",names(state))[1]]] <- E0
            istart <- E0
        } else {
            ## distribute 'E0' value based on dominant eigenvector
            ## here E0 is effectively "number NOT susceptible"
            ee <- round(get_evec(params, testify=testify)*E0)
            if (any(is.na(ee))) {  state[] <- NA; return(state) }
            if (all(ee==0)) {
                if (testify) stop("this case isn't handled for testify")
                ee[["E"]] <- 1  
                warning('initial values too small for rounding')
            }
            istart <- sum(ee)
            state[names(ee)] <- ee
        }
        ## make sure to conserve N by subtracting starting number infected
        ## *after* rounding etc.
        ## FIXME for testify:  (1) make sure get_evec() actually returns appropriate ratios for S
        ##  class; (2) distribute (N-istart) across the S classes, then round
        if (!testify) {
            state[["S"]] <- N-istart
        } else {
            ## if A = testing rate and B = test-return rate then
            ##  du/dt = -A*u + B*n  [where u is untested, n is negative-waiting]
            ##        = -A*u + B*(1-u)  [assuming we're working with proportions]
            ##     -> -u*(A+B) +B =0 -> u = B/(A+B)
            ## FIXME: get_evec() should work for S!
            ufrac <- with(as.list(params),omega/((testing_intensity*W_asymp)+omega))
            state[c("S_u","S_n")] <- round((N-istart)*c(ufrac,1-ufrac))
        }
    } else {
        if (length(names(x))==0) {
            stop("provided state vector must be named")
        }
        if (length(extra_names <- setdiff(names(x),state_names))>0) {
            warning("extra state variables (ignored): ",paste(extra_names,collapse=", "))
        }
        state[names(x)] <- x
    }
    untestify_state <- state ## FIXME: what is this for??
    class(state) <- "state_pansim"
    return(state)
}

##' gradient function for ODE runs
##' @param t time vector
##' @param y state vector
##' @param parms parameter vector
##' @param M rate matrix
gradfun <- function(t, y, parms, M) {
    M <- update_ratemat(M, y, parms)
    foi <- update_foi(y, parms, make_betavec(state=y, parms))
    ## compute 
    flows <- sweep(M, y, MARGIN=1, FUN="*")
    g <- colSums(flows)-rowSums(flows)
    return(list(g,foi=foi))
}

##' Run simulation across a range of times
##' @inheritParams do_step
##' @param nt number of steps to take
##' @param ratemat_args additional arguments to pass to \code{make_ratemat}
##' @param step_args additional arguments to pass to \code{do_step}
##' @param M rate matrix
##' @param use_ode solve as differential equation?
##' @param ode_args additional arguments to ode()
##' @importFrom stats rnbinom
##' @examples
##' params <- read_params("ICU1.csv")
##' r1 <- run_sim_range(params)
##' r2 <- run_sim_range(params,use_ode=TRUE)
##' matplot(r1[,"t"],r1[,-1],type="l",lty=1,log="y")
##' matlines(r2[,"t"],r2[,-1],lty=2)
##' @importFrom dplyr left_join
##' @export
run_sim_range <- function(params
        , state=make_state(params[["N"]], params[["E0"]])
	, nt=100
        , dt=1
        , M = NULL
	, ratemat_args=NULL
        , step_args=NULL
        , use_ode=FALSE
        , ode_args = list()
          ) {
    ## cat("beta0",params[["beta0"]],"\n")
    if (is.null(M)) {
        M <- do.call(make_ratemat,c(list(state=state, params=params), ratemat_args))
    }
    if (use_ode) {
        res <- do.call(deSolve::ode,
                       c(nlist(y=state
                             , times=seq(nt)*dt
                             , func=gradfun
                             , parms = params
                             , M
                               )
                         , ode_args))
        res <- dfs(res)
        if(nrow(res) < nt){
            time_df <- data.frame(time = 1:nt)
            res <- dplyr::left_join(time_df,res)
        }
        names(res)[1] <- "t" ## ode() uses "time"
    } else {
        ## set up output
        foi <- rep(NA,nt)
        res <- matrix(NA, nrow=nt, ncol=length(colnames(M)),
                      dimnames=list(time=seq(nt),
                                    state=colnames(M)))
        ## initialization
        res[1,names(state)] <- state
        if (!has_age(params)) {
            ## FIXME: coherent strategy for accumulating incidence, etc etc
            foi[[1]] <- update_foi(state,params, make_betavec(state, params))
        }
        ## loop
        if (nt>1) {
            for (i in 2:nt) {
                state <- do.call(do_step,
                                 c(nlist(state
                                       , params
                                       , ratemat = M
                                       , dt
                                         )
                                 , step_args))
                if (!has_age(params)) foi[[i]] <- update_foi(state, params, make_betavec(state, params))
                if (!identical(colnames(res),names(state))) browser()
                res[i,] <- state
            }
        }
        res <- dfs(t=seq(nt),res,foi)
    }
    ## need to know true state - for cases with obs error
    attr(res,"state") <- state
    return(res)
}

##' construct a Gamma-distributed delay kernel
##' @param prop area under the curve (proportion reported)
##' @param delay_mean mean value
##' @param delay_cv coeff of var
##' @param max_len maximum delay
##' @param tail_crit criterion for selecting maximum delay
##' @importFrom stats pgamma qgamma
## mean = a*s
## sd = sqrt(a)*s
## cv = 1/sqrt(a)
## s = mean/cv^2
## a = 1/cv^2
make_delay_kernel <- function(prop, delay_mean, delay_cv, max_len=ceiling(tail_val), tail_crit=0.95) {
    
    gamma_shape <- 1/delay_cv^2
    gamma_scale <- delay_mean/gamma_shape
    tail_val <- qgamma(tail_crit, shape=gamma_shape, scale=gamma_scale)
    if (max_len < tail_val) {
        warning(sprintf("max_len (%d) is less than qgamma(%f, %1.1f, %1.1f)=%1.1f",
                        max_len, tail_crit, gamma_shape, gamma_scale, tail_val))
    }
    pp <- diff(pgamma(seq(max_len+1),shape=gamma_shape, scale=gamma_scale))
    pp <- pp/sum(pp) ## normalize to 1
    v <- prop*pp
    return(v)
}

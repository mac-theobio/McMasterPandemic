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
    S <- E <- Ia <- Ip <- Im <- Is <- H  <- NULL
    H2 <- ICUs <- ICUd <- D <- R <- beta0 <- Ca <- Cp  <- NULL
    Cm <- Cs <- alpha <- sigma <- gamma_a <- gamma_m <- gamma_s <- gamma_p  <- NULL
    rho <- delta <- mu <- N <- E0 <- iso_m <- iso_s <- phi1  <- NULL
    phi2 <- psi1 <- psi2 <- psi3 <- c_prop <- c_delaymean <- c_delayCV  <- NULL
    ####
    if (is.null(state)) {
        state <- make_state(N=params[["N"]],E0=1e-3)
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
##'   \item flow diagram: see \url{https://covid-measures.github.io/} 'model details' tab
##'         or \code{../pix/model_schematic.png}
##'   \item parameter definitions: see \code{params_CI_base.csv}, \code{params_ICU_diffs.csv}
##' }
##' 
##' @param state named vector of states
##' @param params named vector of parameters
##' @param do_ICU include additional health utilization compartments
##' @return matrix of (daily) transition rates
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params[["N"]],E0=params[["E0"]])
##' if (require(Matrix)) {
##'   image(Matrix(make_ratemat(state,params)))
##' }
##' @export
make_ratemat <- function(state, params, do_ICU=TRUE) {
    ## circumvent test code analyzers ... problematic ...
    S <- E <- Ia <- Ip <- Im <- Is <- H  <- hosp <- NULL
    H2 <- ICUs <- ICUd <- D <- R <- beta0 <- Ca <- Cp  <- NULL
    Cm <- Cs <- alpha <- sigma <- gamma_a <- gamma_m <- gamma_s <- gamma_p  <- NULL
    rho <- delta <- mu <- N <- E0 <- iso_m <- iso_s <- phi1  <- NULL
    phi2 <- psi1 <- psi2 <- psi3 <- c_prop <- c_delaymean <- c_delayCV  <- NULL
    ## default values, will be masked (on purpose) by unpacking params/state
    nonhosp_mort <- 0
    ####
    np <- length(params)
    ns <- length(state)
    ## make state and param names locally available (similar to with())
    P <- c(as.list(state),as.list(params))
    unpack(P)
    ## blank matrix
    M <- matrix(0,
                nrow=ns, ncol=ns,
                dimnames=list(from=names(state),to=names(state)))
    ## fill entries
    ## NB meaning of iso_* has switched from Stanford model
    Ivec <- c(Ia, Ip, Im,Is)
    Iwt <- beta0/N*c(Ca,Cp,(1-iso_m)*Cm,(1-iso_s)*Cs)
    M["S","E"]   <- sum(Iwt*Ivec)
    M["E","Ia"]  <- alpha*sigma
    M["E","Ip"]  <- (1-alpha)*sigma
    M["Ia","R"]  <- gamma_a
    M["Ip","Im"] <- mu*gamma_p
    M["Ip","Is"] <- (1-mu)*gamma_p
    M["Im","R"]  <- gamma_m
    if (!do_ICU) {
        ## simple hospital model as in Stanford/CEID
        M["Is","H"]  <- gamma_s
        M["H","D"]   <- delta*rho
        M["H","R"]   <- (1-delta)*rho
    } else {
        ## FIXME: A better term than "acute" to mean the opposite of intensive?
        ## four-way split (direct to D, acute care, ICU/survive, ICUD/die)
        M["Is","H"] <- (1-nonhosp_mort)*phi1*gamma_s
        M["Is","ICUs"] <- (1-nonhosp_mort)*(1-phi1)*(1-phi2)*gamma_s
        M["Is","ICUd"] <- (1-nonhosp_mort)*(1-phi1)*phi2*gamma_s
        M["Is","D"] <- nonhosp_mort*gamma_s
        M["ICUs","H2"] <- psi1 ## ICU to post-ICU acute care
        M["ICUd","D"] <- psi2  ## ICU to death
        M["H2","R"]   <- psi3  ## post-ICU to discharge
        ## H now means 'acute care' only; all H survive & are discharged
        M["H","D"]   <- 0
        M["H","R"] <- rho ## all acute-care survive
        if ("hosp" %in% names(state)) {
            M["Is","hosp"] <- M["Is","H"]+M["Is","ICUs"]+M["Is","ICUd"]
            M["hosp","X"] <- 1
            ## assuming that hosp admissions mean *all* (acute-care + ICU)
        }
    }            
    return(M)
}

##' calculate only updated force of infection
##'  at present, this is the only state-dependent \emph{per capita} rate
##'  maybe more efficient than modifying & returning the whole matrix
##' @inheritParams make_ratemat
## FIXME DRY from make_ratemat
update_foi <- function(state, params) {
    ## update infection rate
    with(c(as.list(state),as.list(params)),
         beta0/N*(Ca*Ia+Cp*Ip+(1-iso_m)*Cm*Im+(1-iso_s)*Cs*Is))
}

##' Take a single simulation time step
##' @inheritParams make_ratemat
##' @param ratemat transition matrix
##' @param dt time step (days)
##' @param do_hazard use hazard calculation?
##' @param do_exponential prevent outflow of susceptibles, to create a pure-exponential process?
##' @param stoch_proc stochastic process error?
##' @export
##' @examples
##' params1 <- read_params("ICU1.csv")
##' state1 <- make_state(params=params1)
##' M <- make_ratemat(params=params1, state=state1)
##' s1A <- do_step(state1,params1, M, stoch_proc=TRUE)
do_step <- function(state, params, ratemat, dt=1,
                    do_hazard=FALSE, stoch_proc=FALSE,
                    do_exponential=FALSE) {
    ## FIXME: check (here or elsewhere) for non-integer state and process stoch?
    ## cat("do_step beta0",params[["beta0"]],"\n")
    ratemat["S","E"] <- update_foi(state,params)
    if (!stoch_proc || (!is.null(s <- params[["proc_disp"]]) && s<0)) {
        if (!do_hazard) {
            ## from per capita rates to absolute changes
            flows <- sweep(ratemat, state, MARGIN=1, FUN="*")*dt
        } else {
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
            if (!is.na(proc_disp <- params["proc_disp"])) {
                dW <- pomp::rgammawn(sigma=proc_disp,dt=dt)
            }
            flows[i,-i] <- pomp::reulermultinom(n=1,
                                 size=state[[i]],
                                 rate=ratemat[i,-i],
                                 dt=dW)
            
        }
    }
    outflow <- rowSums(flows)
    if (do_exponential) outflow[["S"]] <- 0
    inflow <-  colSums(flows)
    state <- state - outflow + inflow
    return(state)
}

##' Run pandemic simulation
##' @inheritParams do_step
##' @param start_date starting date (Date or character, any sensible D-M-Y format)
##' @param end_date ending date (ditto)
##' @param params_timevar three-column data frame containing columns 'Date'; 'Symbol' (parameter name/symbol); 'Relative_value' (value \emph{relative to baseline})
##' @param dt time step for do_step
##' @param ratemat_args additional arguments to pass to \code{make_ratemat}
##' @param step_args additional arguments to pass to \code{do_step}
##' @param ndt number of internal time steps per time step
##' @param stoch a logical vector with elements "obs" (add obs error?) and "proc" (add process noise?)
##' @param stoch_start dates on which to enable stochasticity (vector of dates with names 'proc' and 'obs')
##' @param condense condense results?
##' @param condense_args arguments to pass to \code{\link{condense}} (before adding observation error)
##' @examples
##' params <- read_params("ICU1.csv")
##' paramsS <- update(params,c(proc_disp=0.1,obs_disp=100))
##' state <- make_state(params=params)
##' time_pars <- data.frame(Date=c("2020-Mar-20","2020-Mar-25"),
##'                        Symbol=c("beta0","beta0"),
##'                        Relative_value=c(0.7,0.1),
##'                        stringsAsFactors=FALSE)
##' res1 <- run_sim(params,state,start_date="2020-Feb-1",end_date="2020-Jun-1")
##' res1X <- run_sim(params,state,start_date="2020-Feb-1",end_date="2020-Jun-1",condense_args=list(keep_all=TRUE))
##' res1_S <- update(res1, params=paramsS, stoch=c(obs=TRUE, proc=TRUE))
##' res1_t <- update(res1, params_timevar=time_pars)
##' res1_S_t <- update(res1_S, params_timevar=time_pars)
##' res2_S_t <- update(res1_S_t,params=update(paramsS, proc_disp=0.5))
##' res3_S_t <- update(res2_S_t,stoch_start="2020-Apr-1")
##' @importFrom stats rnbinom na.exclude napredict
##' @importFrom anytime anydate
##' @param verbose print messages (e.g. about time-varying parameters)?
##' @export
## FIXME: params_timevar
##   change param name to something less clunky? case-insensitive/partial-match columns? allow Value and Relative_value? (translate to one or the other at R code level, for future low-level code?)
## FIXME: automate state construction better
run_sim <- function(params
        , state=make_state(params[["N"]], params[["E0"]])
        , start_date="2020-Mar-20"
        , end_date="2020-May-1"
        , params_timevar=NULL
        , dt=1, ndt=1  ## FIXME: change default after testing?
        , stoch=c(obs=FALSE,proc=FALSE)
        , stoch_start=NULL
        , ratemat_args=NULL
        , step_args=list()
        , condense = TRUE
        , condense_args=NULL
        , verbose = FALSE
) {
    call <- match.call()

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
    state0 <- state
    nt <- length(date_vec)
    step_args <- c(step_args, list(stoch_proc=stoch[["proc"]]))
    drop_last <- function(x) { x[seq(nrow(x)-1),] }
    M <- do.call(make_ratemat,c(list(state=state, params=params), ratemat_args))
    params0 <- params ## save baseline (time-0) values
    ## no explicit switches, and (no process error) or (process error for full time)
    if (is.null(params_timevar) && (!stoch[["proc"]] || is.null(stoch_start))) {
        switch_times <- NULL
    } else {
        if (!is.null(params_timevar)) {
            ## check column names
            ## FIXME:: tolower()?
            stopifnot(all(c("Date","Symbol","Relative_value") %in%
                          names(params_timevar)))
            params_timevar$Date <- anydate(params_timevar$Date)
        } else {
            params_timevar <- dfs(Date=as.Date(character(0)),Symbol=character(0),Relative_value=numeric(0))
        }
        if (stoch[["proc"]] && !is.null(stoch_start)) {
            params_timevar <- rbind(params_timevar,
                                    dfs(Date=stoch_start[["proc"]],
                                               Symbol="proc_disp",Relative_value=1))
            params[["proc_disp"]] <- -1 ## special value: signal no proc error
        }
        ## convert char to date
        switch_dates <- params_timevar[["Date"]]
        ## match specified times with time sequence
        switch_times <- match(switch_dates, date_vec)
        if (any(is.na(switch_times))) stop("non-matching dates in params_timevar")
    } ## steps

    if (is.null(switch_times)) {
        res <- thin(ndt=ndt,
                    do.call(run_sim_range,
                            nlist(params,state,
                                  nt=nt*ndt,
                                  dt=dt/ndt,M,
                                  ratemat_args,step_args
                                  )
                            ))
    } else {
        t_cur <- 1
        ## want to *include* end date 
        switch_times <- switch_times + 1
        times <- c(1,switch_times,nt+1)
        resList <- list()
        for (i in seq(length(times)-1)) {
            for (j in which(switch_times==times[i])) {
                s <- params_timevar[j,"Symbol"]
                v <- params_timevar[j,"Relative_value"]
                params[[s]] <- params0[[s]]*v
                if (s=="proc_disp") {
                    state <- round(state)
                }
                if (verbose) cat(sprintf("changing value of %s from original %f to %f at time step %d\n",
                            s,params0[[s]],params[[s]],i))
                ## FIXME: so far still assuming that params only change foi
            }
            M["S","E"] <- update_foi(state,params) ## unnecessary?
            resList[[i]] <- drop_last(
                thin(ndt=ndt,
                     do.call(run_sim_range,
                        nlist(params,
                              state,
                              nt=(times[i+1]-times[i]+1)*ndt,
                              dt=dt/ndt,M,
                              ratemat_args,step_args)))
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
    res <- res[,names(res)!="t"]  ## we never want the internal time vector ...
    ## condense here
    if (condense) {
        res <- do.call(condense.pansim,c(list(res,params=params0),condense_args))
    }
    if (stoch[["obs"]]) {
        ## do observation error here
        ## FIXME: allow per-variable obs dispersion; switch to a column-wise operation??
        ## FIXME: warn on mu<0 ? (annoying with ESS machinery)
        m <- res[,-1]
        if (!is.null(stoch_start)) {
            m <- m[res$date>stoch_start[["obs"]],]
        }
        m_rows <- nrow(m)
        mu <- na.exclude(unlist(m))
        mu_S <- rnbinom(length(mu), mu=mu, size=params[["obs_disp"]])
        mu_S <- napredict(attr(mu,"na.action"),mu_S)  ## restore NA values
        res[seq(nrow(res)-m_rows+1,nrow(res)),-1] <- mu_S
    }
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
                       use_eigvec=!is.null(params),
                       params=NULL,
                       x=NULL) {
    ## select vector of state names
    state_names <- switch(type,
                          ## hosp is a hospital-admissions compartment; "X" is a junk compartment
                          ICU1h = c("S","E","Ia","Ip","Im","Is","H","H2","hosp","ICUs","ICUd", "D","R","X"),
                          ICU1 = c("S","E","Ia","Ip","Im","Is","H","H2","ICUs","ICUd", "D","R"),
                          CI =   c("S","E","Ia","Ip","Im","Is","H","D","R"),
                          stop("unknown type")
                          )
    state <- setNames(numeric(length(state_names)),state_names)
    if (is.null(x)) {
        state[["S"]] <- round(N-E0)
        if (!use_eigvec) {
            state[["E"]] <- E0
        } else {
            ## distribute 'E0' value based on dominant eigenvector
            ee <- round(get_evec(params)*E0)
            if (any(is.na(ee))) {  state[] <- NA; return(state) }
            if (all(ee==0)) {
                ee[["E"]] <- 1
                warning('initial values too small for rounding')
            }
            state[names(ee)] <- ee
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
    class(state) <- "state_pansim"
    return(state)
}

##' Run simulation across a range of times
##' @inheritParams do_step
##' @param nt number of steps to take
##' @param ratemat_args additional arguments to pass to \code{make_ratemat}
##' @param step_args additional arguments to pass to \code{do_step}
##' @param M rate matrix
##' @importFrom stats rnbinom
##' @examples
##' params <- read_params("ICU1.csv")
##' run_sim_range(params)
##' @export
run_sim_range <- function(params
        , state=make_state(params[["N"]], params[["E0"]])
	, nt=100
        , dt=1
        , M = NULL
	, ratemat_args=NULL
        , step_args=NULL
          ) {
    ## cat("beta0",params[["beta0"]],"\n")
    if (is.null(M)) {
        M <- do.call(make_ratemat,c(list(state=state, params=params), ratemat_args))
    }
    ## set up output
    foi <- rep(NA,nt)
    res <- matrix(NA, nrow=nt, ncol=length(state),
                  dimnames=list(time=seq(nt),
                                state=names(state)))
    ## initialization
    res[1,] <- state
    foi[[1]] <- update_foi(state,params)
    ## loop
    for (i in 2:nt) {
        state <- do.call(do_step,
                         c(nlist(state
                               , params
                               , ratemat = M
                               , dt
                                 )
                           , step_args))
        foi[[i]] <- update_foi(state, params)
        if (!identical(colnames(res),names(state))) browser()
        res[i,] <- state
    }
    res <- dfs(t=seq(nt),res,foi)
    ## need to know true state - for cases with obs error
    attr(res,"state") <- state
    return(res)
}

##' construct a Gamma-distributed delay kernel
##' @param prop area under the curve (proportion reported)
##' @param delay_mean mean value
##' @param delay_cv coeff of var
##' @param max_len maximum delay
##' @importFrom stats pgamma
## mean = a*s
## sd = sqrt(a)*s
## cv = 1/sqrt(a)
## s = mean/cv^2
## a = 1/cv^2
make_delay_kernel <- function(prop, delay_mean, delay_cv, max_len=10) {
    gamma_shape <- 1/delay_cv^2
    gamma_scale <- delay_mean/gamma_shape
    v <- prop*diff(pgamma(seq(max_len+1),shape=gamma_shape, scale=gamma_scale))
    return(v)
}

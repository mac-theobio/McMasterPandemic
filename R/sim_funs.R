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
    M["S","S"] <- -sum(Ivec*Iwt)
    M["S","Ia"] <- -S*Iwt[["Ia"]]
    M["S","Ip"] <- -S*Iwt[["Ip"]]
    M["S","Im"] <- -S*Iwt[["Im"]]
    M["S","Is"] <- -S*Iwt[["Is"]]
    M["E","Ia"] <- S*Iwt[["Ia"]]
    M["E","Ip"] <- S*Iwt[["Ip"]]
    M["E","Im"] <- S*Iwt[["Im"]]
    M["E","Is"] <- S*Iwt[["Is"]]
    M["E","S"] <- +sum(Ivec*Iwt)
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
    M["R","Ip"] <- gamma_p
    M["R","Im"] <- gamma_m
    M["R","Is"] <- gamma_s
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
    S <- E <- Ia <- Ip <- Im <- Is <- H  <- NULL
    H2 <- ICUs <- ICUd <- D <- R <- beta0 <- Ca <- Cp  <- NULL
    Cm <- Cs <- alpha <- sigma <- gamma_a <- gamma_m <- gamma_s <- gamma_p  <- NULL
    rho <- delta <- mu <- N <- E0 <- iso_m <- iso_s <- phi1  <- NULL
    phi2 <- psi1 <- psi2 <- psi3 <- c_prop <- c_delaymean <- c_delayCV  <- NULL
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
        ## * three-way split (acute care, ICU/survive, ICUD/die)
        M["Is","H"] <- phi1*gamma_s
        M["Is","ICUs"] <- (1-phi1)*(1-phi2)*gamma_s
        M["Is","ICUd"] <- (1-phi1)*phi2*gamma_s
        M["ICUs","H2"] <- psi1 ## ICU to post-ICU acute care
        M["ICUd","D"] <- psi2  ## ICU to death
        M["H2","R"]   <- psi3  ## post-ICU to discharge
        ## H now means 'acute care' only; all H survive & are discharged
        M["H","D"]   <- 0
        M["H","R"] <- rho ## all acute-care survive
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
##' @param stoch stochastic simulation? logical vector for observation and process noise
## (if we do the hazard calculations we can plug in pomp:::reulermultinom()
##  [or overdispersed analogue] directly)
do_step <- function(state, params, ratemat, dt=1,
                    do_hazard=FALSE, stoch=c(obs=FALSE,proc=FALSE),
                    do_exponential=FALSE) {

    ## cat("do_step beta0",params[["beta0"]],"\n")
    ratemat["S","E"] <- update_foi(state,params)
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
    if (!stoch[["proc"]]) {
        outflow <- rowSums(flows)
        if (do_exponential) outflow[["S"]] <- 0
        inflow <-  colSums(flows)
        state <- state - outflow + inflow
    }
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
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params=params)
##' sdate <- "10-Feb-2020" ## arbitrary!
##' time_pars <- data.frame(Date=c("20-Mar-2020","25-Mar-2020"),
##'                        Symbol=c("beta0","beta0"),
##'                        Relative_value=c(0.7,0.1))
##' res <- run_sim(params,state,start_date=sdate,end_date="1-Jun-2020",
##'                    params_timevar=time_pars)
##' summary(res)
##' @importFrom stats rnbinom
##' @param verbose print messages (e.g. about time-varying parameters)?
##' @export
## FIXME: params_timevar
##   change param name to something less clunky? case-insensitive/partial-match columns? allow Value and Relative_value? (translate to one or the other at R code level, for future low-level code?)
## FIXME: automate state construction better
run_sim <- function(params
        , state=make_state(params[["N"]], params[["E0"]])
        , start_date="20-Mar-2020"
        , end_date="1-May-2020"
        , params_timevar=NULL
        , dt=1, ndt=1  ## FIXME: change default after testing?
        , stoch=c(obs=FALSE,proc=FALSE)
        , ratemat_args=NULL
        , step_args=NULL
        , verbose = FALSE
) {
    call <- match.call()

    ## FIXME: *_args approach (specifying arguments to pass through to
    ##  make_ratemat() and do_step) avoids cluttering the argument
    ##  list, but may be harder to translate to lower-level code
    if (dt!=1) warning("nothing has been tested with dt!=1")
    start_date <- ldmy(start_date); end_date <- ldmy(end_date)
    date_vec <- seq(start_date,end_date,by=dt)
    state0 <- state
    nt <- length(date_vec)
    drop_last <- function(x) { x[seq(nrow(x)-1),] }
    M <- do.call(make_ratemat,c(list(state=state, params=params), ratemat_args))
    params0 <- params ## save baseline (time-0) values
    if (is.null(params_timevar)) {
        switch_times <- NULL
    } else {
        ## check column names
        ## FIXME:: tolower()?
        stopifnot(all(c("Date","Symbol","Relative_value") %in%
                      names(params_timevar)))
        ## convert char to date
        switch_dates <- params_timevar[["Date"]] <- ldmy(params_timevar[["Date"]])
        ## match specified times with time sequence
        switch_times <- match(switch_dates, date_vec)
        if (any(is.na(switch_times))) stop("non-matching dates in params_timevar")
    }
    if (is.null(switch_times)) {
        res <- thin(ndt=ndt,
                    do.call(run_sim_range,
                            nlist(params,state,
                                  nt=nt*ndt,
                                  dt=dt/ndt,M,stoch,
                                  ratemat_args,step_args
                                  )
                            ))
    } else {
        t_cur <- 1
        ## want to *include* end date
        times <- c(1,switch_times,nt+1)
        resList <- list()
        for (i in seq(length(times)-1)) {
            for (j in which(switch_times==times[i])) {
                s <- params_timevar[j,"Symbol"]
                v <- params_timevar[j,"Relative_value"]
                params[[s]] <- params0[[s]]*v
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
                              dt=dt/ndt,M,stoch,
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
    res <- data.frame(date=seq(start_date,end_date,by=dt),res)
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
##'     to zero
##' @note \code{"CI"} refers to the Stanford group's
##'     "covid intervention" model.
##' @export
##'
##  FIXME: can pass x, have a name check, fill in zero values
make_state <- function(N=params[["N"]],
                       E0=params[["E0"]],
                       type="ICU1",
                       state_names=NULL,
                       params,
                       x=NULL) {
    ## select vector of state names
    state_names <- switch(type,
       ICU1 = c("S","E","Ia","Ip","Im","Is","H","H2","ICUs","ICUd", "D","R"),
       CI =   c("S","E","Ia","Ip","Im","Is","H","D","R"),
       stop("unknown type")
       )
    state <- setNames(numeric(length(state_names)),state_names)
    if (is.null(x)) {
        state[["S"]] <- N-E0
        state[["E"]] <- E0
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

##' run simulation with specified parameters; extract results
##' matching dates and variable order in data
## FIXME: can this be made into a predict method??
##   with newdata, newparams, newinit
## not sure where to put these ...
##' @param beta0 baseline transmission
##' @param E0 starting value
##' @param data data (for subsetting/matching)
##' @param params parameters
##' @param start_date start date
##' @param values_only return a vector rather than a data frame
##' @importFrom dplyr mutate mutate_at %>% as_tibble right_join
##' @export
predfun <- function(beta0,E0,data,
                    params,
                    start_date="10-Feb-2020",
                    values_only=TRUE) {
    ## global variables
    beta0 <- E0 <- data <- params <- start_date <- values_only <- NULL
    var <- value <- NULL
    ## substitute values into base parameter vector
    params[["beta0"]] <- beta0
    params[["E0"]] <- E0  ## unnecessary?
    state <- make_state(N=params[["N"]],E0=E0) ## assume type==ICU1 for now
    res <- run_sim(params,state,start_date=start_date,
                   end_date=max(data$date)) ## FIXME: pass args?
    ## browser()
    dcomp <- select(data,var,date) %>% mutate_at("var",as.character)
    res2 <- (aggregate(res)
        %>% as_tibble()
        %>% tidyr::pivot_longer(-date,names_to="var")
        %>% right_join(dcomp,by=c("date","var"))
    )
    if (values_only) return(pull(res2,value))
    return(res2)
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
	, stoch=c(obs=FALSE,proc=FALSE)
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
                               , stoch),
                           step_args))
        foi[[i]] <- update_foi(state, params)
        if (!stoch[["obs"]]) {
            res[i,] <- state
        } else {
            res[i,] <- rnbinom(length(state),
                               mu=state,
                               size=params[["obs_disp"]])
        }
    }
    res <- data.frame(t=seq(nt),res,foi)
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

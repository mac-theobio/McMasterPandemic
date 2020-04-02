## FIXME: consistent naming of 'transmat'/'ratemat'/etc.
## FIXME: should eventually split this into multiple files. For now, it's more convenient to have it as one big lump (if/when it becomes a package we can include the tarball and install from source)

##' construct Jacobian matrix for ICU model
##' @param state state vector (named)
##' @param params parameter vector
##' @export
##' @examples
##' params <- read_params(system.file("params","ICU1.csv",package="McMasterPandemic"))
##' state <- make_state(params[["N"]],E=params[["E0"]])
##' ## state[c("E","Ia","Ip","Im","Is")] <- 1
##' state[["E"]] <- 1
##' J <- make_jac(state,params)
##' J["S","S"]
##' Jr <- J[1:6,1:6]
##' round(Jr,3)
##' eigen(Jr)$values
make_jac <- function(state, params) {
    np <- length(params)
    ns <- length(state)
    ## make state and param names locally available (similar to with())
    P <- c(as.list(state),as.list(params))
    attach(P); on.exit(detach(P))
    ## blank matrix
    M <- matrix(0,
                nrow=ns, ncol=ns,
                dimnames=list(from=names(state),to=names(state)))
    Ivec <- c(Ia, Ip, P$Im,Is)
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
    M["E","E"] <- -P$gamma
    M["Ia","E"] <- alpha*P$gamma
    M["Ia","Ia"] <- -lambda_a
    M["Ip","E"] <- (1-alpha)*P$gamma
    M["Ip","Ip"] <- -lambda_p
    M["Im","Ip"] <- mu*lambda_p
    M["Im","Im"] <- -lambda_m
    M["Is","Ip"] <- (1-mu)*lambda_p
    M["Is","Is"] <- -lambda_s
    ## everything else is irrelevant
    return(M)
}

##' Create transition matrix: defines rates (per day) of flow
##' \emph{from} compartment i (row) \emph{to} compartment j (column)
##' base version matches structure of Stanford/Georgia models
##' \itemize{
##' \item flow diagram: see https://covid-measures.github.io/ 'model details' tab or ../pix/model_schematic.png
##' \item parameter definitions: see params_CI_base.csv, params_ICU_diffs.csv
##' }
##' @param state named vector of states
##' @param params named vector of parameters
##' @param do_ICU include additional health utilization compartments
##' @return matrix of (daily) transition rates
##' @examples
##' params <- read_params(system.file("params","ICU1.csv",package="McMasterPandemic"))
##' state <- make_state(params[["N"]],E=params[["E0"]])
##' if (require(Matrix)) {
##'   image(Matrix(make_ratemat(state,params)))
##' }
##' @export
make_ratemat <- function(state, params, do_ICU=TRUE) {
    np <- length(params)
    ns <- length(state)
    ## make state and param names locally available (similar to with())
    P <- c(as.list(state),as.list(params))
    attach(P); on.exit(detach(P))
    ## blank matrix
    M <- matrix(0,
                nrow=ns, ncol=ns,
                dimnames=list(from=names(state),to=names(state)))
    ## fill entries
    ## NB meaning of iso_* has switched from Stanford model
    ## FIXME:: why are gamma(), Im() found rather than the values in P$ ???
    Ivec <- c(Ia, Ip, P$Im,Is)
    Iwt <- beta0/N*c(Ca,Cp,(1-iso_m)*Cm,(1-iso_s)*Cs)
    M["S","E"]   <- sum(Iwt*Ivec)
    M["E","Ia"]  <- alpha*P$gamma
    M["E","Ip"]  <- (1-alpha)*P$gamma
    M["Ia","R"]  <- lambda_a
    M["Ip","Im"] <- mu*lambda_p
    M["Ip","Is"] <- (1-mu)*lambda_p
    M["Im","R"]  <- lambda_m
    if (!do_ICU) {
        ## simple hospital model as in Stanford/CEID
        M["Is","H"]  <- lambda_s
        M["H","D"]   <- delta*rho
        M["H","R"]   <- (1-delta)*rho
    } else {
        ## * three-way split (acute care, ICU/survive, ICUD/die)
        M["Is","H"] <- phi1*lambda_s
        M["Is","ICUs"] <- (1-phi1)*(1-phi2)*lambda_s
        M["Is","ICUd"] <- (1-phi1)*phi2*lambda_s
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
##' @param transmat transition matrix
##' @param dt time step (days)
##' @param do_hazard use hazard calculation?
##' @param stoch stochastic simulation?
## (if we do the hazard calculations we can plug in pomp:::reulermultinom()
##  [or overdispersed analogue] directly)
do_step <- function(state, params, transmat, dt=1,
                           do_hazard=FALSE, stoch=FALSE) {

    transmat["S","E"] <- update_foi(state,params)
    if (!do_hazard) {
        ## from per capita rates to absolute changes
        flows <- sweep(transmat, state, MARGIN=1, FUN="*")*dt
    } else {
        ## use hazard function: assumes exponential change
        ## (constant per capita flows) rather than linear change
        ## (constant absolute flows) within time steps
        ## handle this as in pomp::reulermultinom,
        ## i.e.
        ##    S = sum(r_i)   ## total rate
        ##    p_{ij}=(1-exp(-S*dt))*r_j/S
        ##    p_{ii}= exp(-S*dt)
        S <- rowSums(transmat)
        E <- exp(-S*dt)
        ## prevent division-by-0 (boxes with no outflow) problems (FIXME: DOUBLE-CHECK)
        norm_sum <- ifelse(S==0, 0, state/S)
        flows <- (1-E)*sweep(transmat, norm_sum, MARGIN=1, FUN="*")
        diag(flows) <- 0  ## no flow
    }
    outflow <- rowSums(flows)
    inflow <-  colSums(flows)
    state <- state - outflow + inflow
    return(state)
}

##' Run pandemic simulation
##' @inheritParams do_step
##' @param start_date starting date (Date or character, any sensible D-M-Y format)
##' @param end_date ending date (ditto)
##' @param params_timevar three-column data frame containing columns 'Date'; 'Symbol' (parameter name/symbol); 'Relative_value' (value \emph{relative to baseline})
##' @examples
##' params <- read_params(system.file("params","ICU1.csv",package="McMasterPandemic"))
##' state <- make_state(params[["N"]],E=params[["E0"]])
##' sdate <- "10-Feb-2020" ## arbitrary!
##' time_pars <- data.frame(Date=c("20-Mar-2020","25-Mar-2020"),
##'                        Symbol=c("beta0","beta0"),
##'                        Relative_value=c(0.7,0.1))
##' res <- run_sim(params,state,start_date=sdate,end_date="1-Jun-2020",
##'                    params_timevar=time_pars)
##' summary(res)

##' @export
## FIXME: params_timevar
##   change param name to something less clunky? case-insensitive/partial-match columns? allow Value and Relative_value? (translate to one or the other at R code level, for future low-level code?)
## FIXME: automate state construction better
run_sim <- function(params,
                    state,
                    start_date="20-Mar-2020",
                    end_date="1-May-2020",
                    params_timevar=NULL,
                    dt=1,
                    transmat_args=NULL,
                    step_args=NULL,
                    attach_params=TRUE) {
    call <- match.call()
    ## FIXME: *_args approach (specifying arguments to pass through to
    ##  make_ratemat() and do_step) avoids cluttering the argument
    ##  list, but may be harder to translate to lower-level code
    if (dt!=1) warning("nothing has been tested with dt!=1")
    dfun <- function(x) if (is.character(x)) lubridate::dmy(x) else x
    start_date <- dfun(start_date)
    end_date <- dfun(end_date)
    date_vec <- seq(start_date,end_date,by=dt)
    
    nt <- (as.numeric(end_date-start_date))/dt+1  ## count first date as day 0 (??? FIXME/THINK!)
    ## will non-integer dates work??

    M <- do.call(make_ratemat,c(list(state=state, params=params), transmat_args))
    params0 <- params ## save baseline (time-0) values
    if (is.null(params_timevar)) {
        switch_times <- NULL
    } else {
        ## check column names
        ## FIXME:: tolower()?
        stopifnot(all(c("Date","Symbol","Relative_value") %in%
                      names(params_timevar)))
        ## convert char to date
        switch_dates <- lubridate::dmy(params_timevar[["Date"]])
        ## match specified times with time sequence
        switch_times <- match(switch_dates, date_vec)
        if (any(is.na(switch_times))) stop("non-matching dates in params_timevar")
    }
    ## set up output
    res <- matrix(NA, nrow=nt, ncol=length(state),
                  dimnames=list(time=seq(nt),
                                state=names(state)))
    res[1,] <- state
    for (i in 2:nt) {

        ## apply time-varying/control parameters
        if (i %in% switch_times) {
            for (j in which(switch_times==i)) {
                s <- params_timevar[j,"Symbol"]
                v <- params_timevar[j,"Relative_value"]
                params[[s]] <- params0[[s]]*v
                cat(sprintf("changing value of %s from original %f to %f at time step %d\n",
                            s,params0[[s]],params[[s]],i))
            }
        }
        res[i,] <- do.call(do_step,
                           c(list(state=res[i-1,],
                                  params=params, transmat = M,
                                  dt = dt),
                             step_args))
    }
    res <- data.frame(date=seq(start_date,end_date,by=dt),res)
    ## store everything as attributes
    attr(res,"params") <- params0
    attr(res,"state0") <- state
    attr(res,"start_date") <- start_date
    attr(res,"end_date") <- end_date
    attr(res,"call") <- call
    class(res) <- c("pansim","data.frame")
    return(res)
}

##' @export
print.pansim <- function(x,all=FALSE) {
    if (all) return(unclass(x))
    attr(x,"params") <- NULL
    print(x)
}

    
## packages used 
pkgs <- c("cowplot","tidyverse","ggplot2",
          "igraph", ## for flow chart
          "pomp",   ## for reulermultinom
          "bbmle",
          "emdbook" ## for lambertW
          )

## utility for installing packages
install_pkgs <- function() {
    i1 <- installed.packages()
    pkgs <- pkgs[!pkgs %in% rownames(i1)]
    install.packages(pkgs)
    ## load them all
    suppressMessages(sapply(pkgs,library,character.only=TRUE))
}


## dictionary; internal name, graph label
label_dict <- read.csv(stringsAsFactors=TRUE,
text="
Symbol,Label
S,Susceptible
E,Exposed
I,Total infectious
Ia,Infectious/asymptomatic
Im,Infectious/mild
Is,Infectious/severe
H,Hospital
ICU,ICU
D,Deaths
R,Recovered
")

##' plot method for simulations
##' @param x fitted \code{pansim} object
##' @param drop_states states to \emph{exclude} from plot
##' @param keep_states states to \emph{include} in plot (overrides \code{drop_states})
##' @param aggregate collapse states (e.g. all ICU states -> "ICU") before plotting?  See \code{\link{aggregate.pansim}}
##' @param log plot y-axis on log scale?
##' @export
plot.pansim <- function(x, drop_states=c("S","R","E","I"),
                        keep_states=NULL, aggregate=TRUE,
                        log=FALSE, ...) {
    ## FIXME: check if already aggregated!
    if (aggregate) x <- aggregate(x)
    x <- as_tibble(x)  ## FIXME:: do this upstream?
    if (!is.null(keep_states)) {
        drop_states <- setdiff(names(x), c(keep_states,"date"))
    }
    ## don't try to drop columns that aren't there
    drop_states <- intersect(drop_states,names(x))
    xL <- (x
        %>% as_tibble()
        %>% select(-one_of(drop_states))
        %>% tidyr::pivot_longer(names_to="var", -date)
        %>% mutate(var=forcats::fct_inorder(factor(var)))
    )
    if (log) xL <- dplyr::filter(xL,value>1)
    gg0 <- (ggplot(xL,aes(date,value,colour=var))
        + geom_line()
    )
    if (log) gg0 <- gg0 + scale_y_log10()
    return(gg0)
}


##' retrieve parameters from a CSV file
##' @param fn file name (CSV file containing at least value and symbol columns
##' @param value_col name of column containing values
##' @param symbol_col name of column containing symbols
##' @export
read_params <- function(fn,value_col="Value",symbol_col="Symbol") {
    x <- read.csv(fn,
                  colClasses="character",
                  stringsAsFactors=FALSE,
                  comment="#",
                  na.strings="variable")
    ## evaluate to allow expressions like "1/7" -> numeric
    x[[value_col]] <- vapply(x[[value_col]], function(z) eval(parse(text=z)), numeric(1))
    res <- setNames(x[[value_col]],x[[symbol_col]])
    class(res) <- "panparams"
    return(res)
}

##' write parameters to CSV file
##' @param fn file name
##' @param params a params object
##' @param label a label for the parameters
##' @export
write_params <- function(params, fn, label) {
    writeLines(con=fn,
           c(paste("#",label),
             sprintf("# Date: %s",format(Sys.time(),"%d %b %Y"))))
    ## unavoidable warning "appending column names to file"
    suppressWarnings(write.table(data.frame(
        Symbol=names(params),
        Value=unclass(params)), file=fn,
        row.names=FALSE, append=TRUE,
        sep=","))
}

##' Collapse columns (infected, ICU, hospitalized) in a pansim output
##' @param x a pansim object
##' @export
aggregate.pansim <- function(x) {
    ## FIXME: less clunky way to do this? Collapse columns *if* present
    ##   but account for the fact that they might be missing in some variants
    ## FIXME: extend to aggregate, S, E, etc. as we add space / testing / age
    ## may need to go tidy?
    c0 <- class(x)
    ## collapse columns and add, if present
    add_col <- function(dd,name,regex) {
        vars <- grep(regex, names(x), value=TRUE)
        if (length(vars)>0) {
            dd[[name]] <- rowSums(x[vars])
        }
        return(dd)
    }
    dd <- x[c("date","S","E")]
    dd <- add_col(dd,"I","^I[^C]")
    dd <- add_col(dd,"H","^H")
    dd <- add_col(dd,"ICU","^ICU")
    dd <- data.frame(dd,R=x[["R"]])
    dd <- add_col(dd,"discharge","discharge")
    dd <- data.frame(dd,D=x[["D"]])
    class(dd) <- c0 ## make sure class is restored
    return(dd)
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

get_r <- function(params) {
    ## will have to be done by constructing Jacobian?
} 


##' @export
summary.pansim <- function(x, ...) {
    ## FIXME: get ventilators by multiplying ICU by 0.86?
    ## FIXME: prettier?
    xa <- aggregate(x)
    attach(xa); on.exit(detach(xa))
    res <- data.frame(peak_ICU_date=date[which.max(ICU)],
             peak_ICU_val=round(max(ICU)),
             peak_H_date=date[which.max(H)],
             peak_H_val=round(max(H)))
    ## FIXME: report time-varying R0
    if (!is.null(p <- attr(x,"params"))) {
        res <- data.frame(res,R0=get_R0(p))
    }
    class(res) <- c("summary.pansim","data.frame")
    res
}

##' generate initial state vector
##' @param N population size
##' @param E0 initial number exposed
##' @param type (character) specify what model type this is intended for; determines state names
##' @param state_names vector of state names, must include S and E
##' @export
make_state <- function(N,E0,type="ICU1",state_names=NULL) {
    state_names <- switch(type,
       ICU1 = c("S","E","Ia","Ip","Im","Is","H","H2","ICUs","ICUd", "D","R"),
       CI =   c("S","E","Ia","Ip","Im","Is","H","D","R"),
       stop("unknown type")
       )
    state <- setNames(numeric(length(state_names)),state_names)
    state[["S"]] <- N-E0
    state[["E"]] <- E0
    class(state) <- "state.pansim"
    return(state)
}

##' run simulation with specified parameters; extract results
##' matching dates and variable order in data
## FIXME: can this be made into a predict method??
##   with newdata, newparams, newinit

##' @export
predfun <- function(beta0,E0,data,
                    params,
                    start_date="10-Feb-2020",
                    values_only=TRUE) {
    require("dplyr")
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


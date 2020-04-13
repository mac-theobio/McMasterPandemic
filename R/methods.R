##' @export
print.pansim <- function(x,all=FALSE,...) {
    ## FIXME: is this the best way?
    ## use tibbles or not?
    class(x) <- "data.frame"
    if (all) print(x)
    attr(x,"params") <- NULL
    print(x)
}

calc_reports <- function(x,params) {
    ## compute incidence and reports (as convolution of incidence)
    incidence <- x$foi*x$S
    unpack(as.list(params))
    kern <- make_delay_kernel(c_prop,
                              c_delay_mean,
                              c_delay_cv,
                              max_len=10)
    ## FIXME: don't hard-code max len ...
    report <- as.numeric(stats::filter(incidence,kern,
                                       sides=1))

    return(data.frame(incidence, report))
}

## FIXME: allow faceting automatically? (each var alone or by groups?)
## don't compare prevalences and incidences?
##' plot method for simulations
##' @param x fitted \code{pansim} object
##' @param drop_states states to \emph{exclude} from plot
##' @param keep_states states to \emph{include} in plot (overrides \code{drop_states})
##' @param aggregate collapse states (e.g. all ICU states -> "ICU") before plotting?  See \code{\link{aggregate.pansim}}
##' @param log plot y-axis on log scale?
##' @param show_times indicate times when parameters changed?
##' @param ... additional arguments to \code{\link{aggregate.pansim}}
##' @importFrom ggplot2 ggplot geom_line aes geom_vline scale_y_log10
##' @importFrom dplyr one_of
##' @export
plot.pansim <- function(x, drop_states=c("t","S","R","E","I","incidence"),
                        keep_states=NULL, aggregate=TRUE,
                        log=FALSE, show_times=TRUE, ...) {
    ## global variables
    var <- value <- NULL
    ## attributes get lost somewhere below ...
    ptv <- attr(x,"params_timevar")
    if (aggregate && !isTRUE(attr(x,"aggregated"))) {
        x <- aggregate(x)
    }
    x <- as_tibble(x)  ## FIXME:: do this upstream?
    if (!is.null(keep_states)) {
        drop_states <- setdiff(names(x), c(keep_states,"date"))
    }
    ## don't try to drop columns that aren't there
    ## FIXME: use aggregate.pansim method?
    drop_states <- intersect(drop_states,names(x))
    ## FIXME: don't pivot if already pivoted
    xL <- (x
        %>% as_tibble()
        %>% dplyr::select(-one_of(drop_states))
        %>% tidyr::pivot_longer(names_to="var", -date)
        %>% mutate(var=forcats::fct_inorder(factor(var)))
    )
    if (log) xL <- dplyr::filter(xL,value>=1)
    gg0 <- (ggplot(xL,aes(date,value,colour=var))
        + geom_line()
    )
    if (log) gg0 <- gg0 + scale_y_log10()
    if (show_times && !is.null(ptv)) {
        gg0 <- gg0 + geom_vline(xintercept=ptv$Date,lty=2)
    }
    return(gg0)
}

##' Collapse columns (infected, ICU, hospitalized) in a pansim output
##' @param x a pansim object
##' @param agg_states aggregate states (and add case reports)?
##' @param pivot return long-format tibble instead of wide data frame?
##' @param keep_vars variables to retain (in addition to date) if pivoting
##' @param t_agg_start starting date for temporal aggregation
##' @param t_agg_period time period for temporal aggregation (e.g. "7 days", see \code{\link{seq.Date}})
##' @param t_agg_fun temporal aggregation function (applied across all variables) \emph{or} a list of the form \code{list(FUN1=c('var1','var2'), FUN2=c('var3', 'var4'))} (temporal aggregation is done after state aggregation, so variable names specified should be adjusted appropriately)
##' @param add_reports add incidence and case reports?
##' @param ... unused, for generic consistency
##' @importFrom stats aggregate
##' @importFrom dplyr %>% as_tibble
##' @importFrom tidyr pivot_longer
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params=params)
##' sdate <- "10-Feb-2020" ## arbitrary!
##' res <- run_sim(params,state,start_date=sdate,end_date="1-Jun-2020")
##' a1 <- aggregate(res, t_agg_start="12-Feb-2020",t_agg_period="7 days",t_agg_fun=sum, agg_state=FALSE)
##' plot(a1) + ggplot2::geom_point()
##' ## column-specific aggregation
##' first <- dplyr::first
##' a2 <- aggregate(res, t_agg_start="12-Feb-2020",t_agg_period="7 days",
##'         t_agg_fun=list(mean=c("H","ICU","I"),
##'                first=c("D"),sum=c("report")))
##' @export
aggregate.pansim <- function(x,pivot=FALSE,keep_vars=c("H","ICU","D","report"),
                             agg_states=TRUE,
                             add_reports=TRUE,
                             t_agg_start=NULL,
                             t_agg_period=NULL,
                             t_agg_fun=mean,
                             ...) {
    ## FIXME: less clunky way to do this? Collapse columns *if* present
    ##   but account for the fact that they might be missing in some variants
    ## FIXME: extend to aggregate, S, E, etc. as we add space / testing / age
    ## may need to go tidy?
    ## global variables
    c_prop <- c_delay_mean <- c_delay_cv <- NULL
    c0 <- class(x)
    dd <- x
    if (agg_states) {
        ## collapse columns and add, if present
        add_col <- function(dd,name,regex) {
            vars <- grep(regex, names(x), value=TRUE)
            if (length(vars)>0) {
                dd[[name]] <- rowSums(x[vars])
            }
            return(dd)
        }
        ## Sorry, Ben; no good way to get reports and Is, Im
        ## Add reports as a separate argument instead of part of agg_states?
        ## MLi: I am sorry too, we needed more states 2:27am.
        dd <- dd[c("date","S","E","Ia", "Is", "Im", "H","H2","ICUs","ICUd","foi")]
        dd <- add_col(dd,"I","^I[^C]")
        dd <- add_col(dd,"H","^H")
        dd <- add_col(dd,"ICU","^ICU")
        dd <- data.frame(dd,R=x[["R"]])
        dd <- add_col(dd,"discharge","discharge")
        dd <- data.frame(dd,D=x[["D"]])
        params <- attr(x,"params")
    } ## if agg_states
    if (add_reports) {
        if (!"c_delay_mean" %in% names(params)) {
            warning("add_reports requested but delay parameters missing")
        } else {
            dd <- data.frame(dd, calc_reports(x,params))
        }
    }
    ## do temporal aggregation, if requested
    if (!is.null(t_agg_start)) {
        agg_datevec <- seq.Date(ldmy(t_agg_start),max(dd$date)+30,by=t_agg_period)
        agg_period <- cut.Date(dd$date,agg_datevec)
        ## set to *last* day of period
        ap <- as.Date(levels(agg_period))
        dt <- as.numeric(diff(ap))[1]
        levels(agg_period) <- as.character(ap+dt-1)
        if (is.function(t_agg_fun)) {
            dd <- stats::aggregate.data.frame(dplyr::select(dd,-date),
                                              by=list(date=agg_period),
                                              FUN=t_agg_fun)
        } else {
            if (!is.list(t_agg_fun)) {
                stop("t_agg_fun should be either a single function or a list of the form list(FUN1=c('var1','var2'), FUN2=c('var3', 'var4'))")
            }
            dd_tmp <- list()
            for (i in seq_along(t_agg_fun)) {
                FUN <- get(names(t_agg_fun)[i])
                for (j in seq_along(t_agg_fun[[i]])) {
                    var <- t_agg_fun[[i]][[j]]
                    dd_tmp <- c(dd_tmp,
                                setNames(list(stats::aggregate(dd[[var]],
                                                               by=list(date=agg_period),
                                                               FUN=FUN)[,2]),var))
                }
            }
            dd <- do.call(data.frame,c(list(date=unique(na.omit(agg_period))),dd_tmp))
            ## FIXME: fix order of columns?
        }
        dd$date <- as.Date(dd$date)
    }
    class(dd) <- c0 ## make sure class is restored
    if (!pivot) return(dd)
    ## OTHERWISE long form: more convenient for regressions etc.
    ## keep_vars <- names(dd)
    keep_vars <- intersect(keep_vars,names(dd))
    dd <- (dd
        %>% as_tibble()
        %>% dplyr::select(c("date",keep_vars))
        %>% pivot_longer(names_to="var",-date)
    )
    ## FIXME: distinguish between agg w/ and w/o pivot?
    attr(dd,"aggregated") <- TRUE
    ## FIXME: restore other attributes!
    class(dd) <- c0 ## restore class
    return(dd)
}

##' summary method for \code{pansim} parameter objects
##' @param object a params vector
##' @param ... unused
##' @export
summary.params_pansim <- function(object, ...) {
    ## FIXME: include kappa once we know what works
    ## (analytical vs kernel vs ...)
    res <- c(r0=get_r(object),R0=get_R0(object),Gbar=get_Gbar(object))
    res["dbl_time"] <- log(2)/res["r0"]
    return(res)
}

##' @export
summary.pansim <- function(object, ...) {
    ## global variables
    ICU <- H <- NULL
    ## FIXME: get ventilators by multiplying ICU by 0.86?
    ## FIXME: prettier?
    xa <- aggregate(object)
    unpack(xa)
    res <- data.frame(peak_ICU_date=xa$date[which.max(ICU)],
             peak_ICU_val=round(max(ICU)),
             peak_H_date=xa$date[which.max(H)],
             peak_H_val=round(max(H)))
    ## FIXME: report time-varying R0
    if (!is.null(p <- attr(object,"params"))) {
        res <- data.frame(res,R0=get_R0(p))
    }
    class(res) <- c("summary.pansim","data.frame")
    res
}

## hard-coded param descriptions
## use as backup if params don't have them attached
param_meanings <- c(
    beta0 = "transmission rate",
    Ca = "relative asymptomatic transmissibility",
    Cp = "relative pre-symptomatic transmissibility",
    Cs = "relative severely symptomatic transmissibility",
    Cm = "relative mildly symptomatic transmissibility",
    alpha = "proportion of infections that are asymptomatic",
    sigma = "1 / mean LATENT period",
    gamma_a = "1 / mean days in asymptomatic infectious class",
    gamma_s = "1 / mean days in severely symptomatic infectious class",
    gamma_m = "1 / mean days in mildly symptomatic infectious class",
    gamma_p = "1 / mean days in pre-symptomatic infectious class",
    rho = "1 / mean days in acute care",
    delta = "proportion of acute care patients who die",
    mu = "proportion of symptomatic infections that are mild",
    N = "total population size",
    E0 = "number of initially exposed individuals",
    iso_m = "proportion of mildly symptomatic patients who are isolated",
    iso_s = "proportion of severely symptomatic patients who are isolated",
    phi1 = "proportion of severe infections that do NOT require ICU",
    phi2 = "proportion of ICU patients who die",
    psi1 = "1 / mean days in ICU if survive",
    psi2 = "1 / mean days in ICU if die",
    psi3 = "1 / mean days post-ICU until discharge"
)

##' Describe parameters
##'
##' Create a data frame with symbols, values and meanings of parameters
##' @param x a \code{params_pansim} object
##' @export
describe_params <- function(x) {
    if (!is.null(attr(x,"description"))) {
        x_meanings <- attr(x,"description")[names(x)]
    } else {  ## backup/built-in
        x_meanings <- param_meanings[names(x)]
    }
    xout <- data.frame(symbol=names(x),
                       ##value=round(as.numeric(x),3),
                       value=sprintf("%.3g", as.numeric(x)),
                       meaning=x_meanings)
    rownames(xout) <- NULL ## redundant
    return(xout)
}

##' print parameters, possibly with detailed description
##'
##' If the detailed description is requested, it is returned as a data frame.
##' @param x an object of class \code{params_pansim} (parameters for pandemic simulation)
##' @param describe print full description?
##' @param ... (unused, for generic consistency)
##' @examples
##' params <- read_params(system.file("params","ICU1.csv",package="McMasterPandemic"))
##' print(params)
##' print(params,describe=TRUE)
##' @export
## FIXME: prettier printing, e.g. detect "1/" or "proportion"
print.params_pansim <- function( x, describe=FALSE, ... ) {
    if (!describe) {
        attr(x,"description") <- NULL
        print(unclass(x))
    } else {
        print(describe_params(x))
    }
}

##' @export
update.params_pansim <- function(object, ...) {
    L <- list(...)
    nm <- names(L)
    for (i in seq_along(L)) {
        object[[nm[i]]] <- L[[i]]
    }
    return(object)
}

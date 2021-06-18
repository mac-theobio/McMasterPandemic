##' @export
print.pansim <- function(x,all=FALSE,...) {
    ## FIXME: is this the best way?
    ## use tibbles or not?
    class(x) <- "data.frame"
   if (all) print(x)
    attr(x,"params") <- NULL
    print(x)
}

##' calculate convolution
##' @param i an incidence time series
##' @param params a list or vector containing elements \code{c_prop}, \code{c_delay_mean}, \code{c_delay_cv}
##' @export
calc_conv <- function(i,params) {
    c_prop <- c_delay_mean <- c_delay_cv <- NULL
    kern <- with(as.list(params),
                 make_delay_kernel(c_prop,
                                   c_delay_mean,
                                   c_delay_cv)
                 )
    ret <- as.data.frame(as.numeric(stats::filter(i,kern,sides=1)))
    ## if parameters are ageified, keep age-stratified reports in output too
    if(has_age(params)){
        state_suffixes <- sub("^incidence", "",
                              grep("^incidence_", names(i), value = TRUE))
        names(ret) <- paste0("report", state_suffixes)
        ## add total reports
        ret$report <- rowSums(ret)
    } else{
        names(ret) <- c("report")
    }
    return(ret)
}

## compute incidence and reports (as convolution of incidence)
calc_reports <- function(x, params, add_cumrep=FALSE) {
    ## FIXME: dt==1 !
    ## calculate incidence (for each age group, if it exists)
    if(length(grep("^S_", names(x), value = TRUE)) > 0){
        state_suffixes <- sub("^S", "", grep("^S_", names(x), value = TRUE))
        incidence <- lapply(state_suffixes,
                            function(suff) x[[paste0("foi", suff)]]*x[[paste0("S", suff)]])
        ## convert to df
        incidence <- as.data.frame(do.call(cbind, incidence))
        names(incidence) <- paste0("incidence", state_suffixes)
    } else {
        ## otherwise, just do total incidence
        incidence <- data.frame(incidence = x$foi*x$S)
    }

    ## add total incidence to output
    incidence$incidence <- rowSums(incidence)

    ## FIXME: only calculates total reports right now, not age-stratified see
    ## what happens if, within calc_conv, we do a convolution for eacha age
    ## group and then sum up vs one convolution over total incidence
    report <- calc_conv(incidence$incidence,params)

    ret <- dfs(incidence, report)
    ## FIXME: take out the cum rep stuff?  this is the wrong place for it,
    ## it's generally happening *before* the addition of obs error
    if (add_cumrep) {
        cumRep <- cumsum(ifelse(!is.na(report), report, 0))
        ret <- data.frame(ret, cumRep)
    }
    return(ret)
}

#' Prepare age-structured simulation results data frame for plotting
#'
#' @param res age-structured simulation result
#' @param split_age_vax should the age and vaccination category be split into different columns?
#' @inheritParams plot.pansim
#' @importFrom forcats as_factor
#' @importFrom stringr str_replace
#'
#' @return
#' @export
prep_res_for_plotting <- function(res,
                                  drop_states = NULL,
                                  condense_I = FALSE,
                                  split_age_vax = FALSE){

    if(has_vax(res) && split_age_vax){
        into <- c("state", "age", "vaccination")
    } else {
        into <- c("state", "age")
    }

    (res
     %>% select(-starts_with("foi"))
     %>% pivot_longer(-date)
     %>% separate(name, into = into,
                  sep = "_", extra = "merge")
    ) -> res

    ## condense I cats
    if(condense_I){
        (res
         ## convert state column to factor to maintain original order of variables
         %>% mutate(state = as_factor(str_replace(state,
                                                           "I[amps]", "I")))
         %>% group_by(across(c(-value)))
         %>% summarise(value = sum(value), .groups = "drop")
        ) -> res
    }

    (res
        ## fix age labels
        %>% mutate(age = str_replace(age, "\\.$", "\\+"))
        %>% mutate(age = str_replace(age, "\\.", "-"))
        ## convert categories to factors to maintain ordering
        %>% mutate(across(where(is.character), ~ as_factor(.)))
    ) -> res

    if(!is.null(drop_states)){
        res <- res %>% filter(!(state %in% drop_states))
    }

    return(res)
}

#' Plot age-structured simulation result faceted by age categories
#'
#' @param res age-structured simulation result
#' @inheritParams plot.pansim
#' @importFrom dplyr vars
#' @importFrom ggplot2 scale_x_date
#'
#' @return
#' @export
plot_res_by_age <- function(res, drop_states = NULL,
                            condense_I = FALSE,
                            show_times = TRUE){
    ## get time-varying params attribute, if it exists
    ptv <- attr(res,"params_timevar")

    (prep_res_for_plotting(res, drop_states, condense_I)
     %>% ggplot(aes(x = date, y = value, colour = state))
     + geom_line()
     + facet_wrap(vars(age))
     + scale_x_date(date_breaks = "1 month",
                    date_labels = "%b")
     # + scale_y_continuous(labels = scales::label_number_si())
    ) -> gg

    if (show_times && !is.null(ptv)) {
        gg <- gg + geom_vline(xintercept=ptv$Date,lty=2)
    }

    return(gg)
}

#' Plot age-structured simulation result faceted by state
#'
#' @param res age-structured simulation result
#' @inheritParams plot.pansim
#'
#' @return
#' @export
plot_res_by_state <- function(res, drop_states = NULL,
                              condense_I = FALSE,
                              show_times = TRUE){
    ## get time-varying params attribute, if it exists
    ptv <- attr(res,"params_timevar")

    if(has_vax(res)){
      plot_setup <- (prep_res_for_plotting(res, drop_states, condense_I, split_age_vax = TRUE) %>%
      ggplot(aes(x = date, y = value, colour = age, linetype = vaccination)))
    } else {
      plot_setup <- (prep_res_for_plotting(res, drop_states, condense_I) %>%
                       ggplot(aes(x = date, y = value, colour = age)))
    }

    (plot_setup
     + geom_line()
     + facet_wrap(vars(state), scales = "free_y")
     + scale_x_date(date_breaks = "1 month",
                    date_labels = "%b")
     # + scale_y_continuous(labels = scales::label_number_si())
    ) -> gg

    if (show_times && !is.null(ptv)) {
        gg <- gg + geom_vline(xintercept=ptv$Date,lty=2)
    }

    return(gg)
}

## FIXME: allow faceting automatically? (each var alone or by groups?) don't
## compare prevalences and incidences?
## FIXME: incorporate age-specific plotting in setup of base plot,
## then add stuff (like vlines for timepars) overtop,
## just in this plotting function (not separately in each plot style) currently,
## this is done in a redundant way currently
##' plot method for simulations
##' @param x fitted \code{pansim} object
##' @param drop_states states to \emph{exclude} from plot
##' @param keep_states states to \emph{include} in plot (overrides \code{drop_states})
##' @param condense condense states (e.g. all ICU states -> "ICU") before plotting?  See \code{\link{condense.pansim}}
##' @param facet_by_age if this is an age-structured simulation, do we want to facet by age? if FALSE, facet by state (default)
##' @param log plot y-axis on log scale?
##' @param log_lwr lower bound for log scale
##' @param show_times indicate times when parameters changed?
##' @param ... additional arguments to \code{\link{condense.pansim}}
##' @importFrom ggplot2 ggplot geom_line aes geom_vline scale_y_log10 geom_ribbon
##' @importFrom dplyr one_of
##' @return a \code{\link[ggplot2]{ggplot}} object
##' @export
plot.pansim <- function(x, drop_states=c("t","S","R","E","I","X","incidence"),
                        keep_states=NULL, condense=FALSE,
                        facet_by_age = FALSE,
                        log=FALSE,
                        log_lwr=1,
                        show_times=TRUE, ...) {
    ## global variables
    var <- value <- NULL

    ## if age-structured, use a different plotting method
    if(has_age(x)){
        plot_args <- list(res = x,
                          drop_states = drop_states,
                          condense_I = condense,
                          show_times = show_times)
        if(facet_by_age){
            return(do.call(plot_res_by_age, plot_args))
        } else{
            return(do.call(plot_res_by_state, plot_args))
        }
    }

    ## attributes get lost somewhere below ...
    ptv <- attr(x,"params_timevar")

    if (!is.null(keep_states)) {
        drop_states <- setdiff(names(x), c(keep_states,"date"))
    }
    ## don't try to drop columns that aren't there
    ## FIXME: use condense.pansim method?
    drop_states <- intersect(drop_states,names(x))
    ## FIXME: don't pivot if already pivoted
    xL <- (pivot(x)
        %>% dplyr::filter(!var %in% drop_states)
        %>% dplyr::mutate(var=factor(var,levels=unique(var)))
        ## FIXME: order factor in pivot?
    )
    if (log) xL <- dplyr::filter(xL,value>=log_lwr)
    gg0 <- (ggplot(xL,aes(date,value,colour=var))
        + geom_line()
    )
    if (log) gg0 <- gg0 + scale_y_log10()
    if (show_times && !is.null(ptv)) {
        gg0 <- gg0 + geom_vline(xintercept=ptv$Date,lty=2)
    }
    return(gg0)
}



##' condense an object
##' @param object an object to condense
##' @param ... additional arguments
##' @export
condense <- function (object, ...)  {
    UseMethod("condense")
}

##' pivot an object
##' @param object an object to pivot
##' @param ... additional arguments
##' @export
pivot <- function (object, ...)  {
    UseMethod("pivot")
}

##' @export
##' @importFrom dplyr %>%
pivot.pansim <- function(object, ...) {
    check_dots(...)
    dd <- (object
        %>% dplyr::as_tibble()
        %>% tidyr::pivot_longer(names_to="var",-date)
    )
    return(dd)
}

## test whether variables have already been condensed
is_condensed <- function(x) "I" %in% names(x)
has_report <- function(x) "report" %in% names(x)

##' Condense columns (infected, ICU, hospitalized) in a pansim output
##' @param object a pansim object
##' @param add_reports add incidence and case reports?
##' @param diff_deaths compute first differences of death series to get daily deaths?
##' @param keep_all keep unaggregated variables in data frame as well?
##' @param cum_reports compute cumulative reports?
##' @param params parameters (for defining convolution kernel for reports, and zeta and N for het-S)
##' @param het_S compute hetS = (S/N)^(1+zeta) ?
##' @param ... additional args
##' @export
condense.pansim <-  function(object, add_reports=TRUE,
                             diff_deaths=TRUE,
                             cum_reports=FALSE,
                             het_S=FALSE,
                             keep_all=FALSE,
                             params = attr(object,"params"),
                             ...)
{
    check_dots(...)
    aa <- get_attr(object)

    regex_expanded_suffix <- case_when(
      has_age(params) & has_vax(params) ~ "_[0-9]+\\.[0-9]?[0-9]?_",
      has_age(params) ~ "_[0-9]+\\.[0-9]?[0-9]?$",
      has_vax(params) ~ "vax",
      T ~ "")

    ## condense columns and add, if present
    add_col <- function(dd, name, regex, collapse = TRUE) {
        vars <- grep(regex, names(object), value=TRUE)
        if (length(vars)>0) {
            if (collapse){
                dd[[name]] <- rowSums(object[vars])
            } else {
                dd <- dfs(dd, object[vars])
            }

        }
        return(dd)
    }

    ## should incidence be added for expanded sims?
    ## (need foi columns, if so)
    add_expanded_reports <- (add_reports
                              && (has_age(params) || has_vax(params))
                              && any(grepl("^foi", names(object))))

    ## check first if condensed
    ## FIXME: rearrange logic?
    if (is_condensed(object)) {
        dd <- object
    } else {
        if (keep_all) {
            dd <- object
        } else {
            ## start by keeping date
            dd <- object["date"]

            ## if foi by subcategory is in results,
            ## keep unaggregated S categories for now
            ## (so that we can calculate incidence)
            if (add_expanded_reports){
                dd <- add_col(dd,"S","^S", collapse = FALSE)
            }

            ## keep susc and exposed classes
            for (n in c("S","E")) {
                dd <- add_col(dd, n, paste0("^",n))
            }

            dd <- add_col(dd,"I","^I[^C]") ## collapse all I* variables that aren't ICU
            for (n in c("H", "ICU","R")) {
                dd <- add_col(dd,n,paste0("^",n))
            }

            diff_vars <- c(X="hosp",D="death",N="negtest",P="postest")
            for (i in seq_along(diff_vars)) {
                nm <- names(diff_vars)[i]
                re <- paste0("^",nm)
                if (any(grepl(re,names(object)))) {
                    tot <- rowSums(object[grepl(re,names(object))])
                    dd[[diff_vars[i]]] <- c(NA,diff(tot))
                    dd <- add_col(dd,nm,re)
                }
            }
            ## keep foi (if it exists) as a single column
            if ("foi" %in% names(object)){
                dd <- data.frame(dd, foi = object[["foi"]])
            } else {
                ## if foi is expanded, keep without condensing
                ## in case we want to calculate incidence later
                if (any(grepl("^foi", names(object)))){
                    dd <- add_col(dd,"foi","^foi", collapse = FALSE)
                }
            }

        }  ## not keep_all
    } ## already condensed

    ## no check if we should add reports
    if (add_reports) {
        ## if reports aren't already in the output object
        if(!has_report(dd)){
            ## and they're not already in the input object
            if(!has_report(object)){
                if (!"c_delay_mean" %in% names(params)) {
                    warning("add_reports requested but delay parameters missing")
                } else {
                    if (cum_reports) warning("cum_reports is deprecated (reports are cumulated in run_sim)")
                    ## if we don't either have condensed S or age-specific S, can't generate reports
                    # if (!any("S" %in% names(dd), grepl("^S_[0-9]+", names(dd)))) {
                    #     stop("need either condensed S or age-specific S to compute reports (at least one is missing)")
                    # }
                    cr <- calc_reports(dd, params, add_cumrep=cum_reports)
                    dd <- data.frame(dd, cr)
                }
            } else {
                ## the reports are in the input object, so just copy them over
                dd <- dfs(dd, object[grepl("^(incidence|report|cumRep)",
                                           names(object))])
            }
        }
        ## but if not keep_all, drop any expanded reports (if they've been
        ## added to the output object)
        if(!keep_all && (has_age(params) || has_vax(params))){
            dd <- dd[!grepl(paste0("^(incidence|report)", regex_expanded_suffix),
                            names(dd))]
        }
    } ## add_reports

        if (het_S) {
            dd$hetS <- (dd$S/params[["N"]])^(1+params[["zeta"]])
        }

        ## drop age-specific S classes
        ## (if they still exist from the incidence calc)
        if(add_expanded_reports && !keep_all){
            dd <- dd[!grepl("^S_", names(dd))]
        }
        dd <- put_attr(dd,aa)
        return(dd)
}

##' Temporal aggregation of pansim objects
##' @param x a pansim object
##' @param start starting date for temporal aggregation
##' @param period time period for temporal aggregation (e.g. "7 days", see \code{\link{seq.Date}})
##' @param FUN temporal aggregation function (applied across all variables) \emph{or} a list of the form \code{list(FUN1=c('var1','var2'), FUN2=c('var3', 'var4'))} (temporal aggregation is done after state aggregation, so variable names specified should be adjusted appropriately)
##' @param extend how far to extend time series for aggregation purposes
##' @param fixed_vars treat names in FUN as fixed strings rather than regular expressions?
##' @param ... unused, for generic consistency
##' @importFrom stats aggregate
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params=params)
##' sdate <- "2020-02-10" ## arbitrary!
##' res <- run_sim(params,state,start_date=sdate,end_date="2020-06-01")
##' a1 <- aggregate(res, start="2020-02-12",period="7 days", FUN=sum)
##' ## column-specific aggregation
##' first <- dplyr::first
##' agg_funs <- list(mean=c("H","ICU","I"), sum=c("report","death"))
##' a2 <- aggregate(condense(res), start="2020-02-12",period="7 days", FUN=agg_funs)
##' @export
aggregate.pansim <- function(x,
                             start=NULL,
                             period=NULL,
                             FUN=mean,
                             fixed_vars=TRUE,
                             extend=30,
                             ...) {
    check_dots(...)
    aa <- get_attr(x)
    dd <- x
    ## start <- agg_list[["t_agg_start"]]
    ## period <- agg_list[["t_agg_period"]]
    ## FUN <- agg_list$t_agg_fun
    agg_datevec <- seq.Date(as.Date(start),max(dd$date)+extend,
                            by=period)
    agg_period <- cut.Date(dd$date,agg_datevec)
    ## set to *last* day of period
    ap <- as.Date(levels(agg_period))
    dt <- as.numeric(diff(ap))[1]
    levels(agg_period) <- as.character(ap+dt-1)
    if (is.function(FUN)) {
        dd <- stats::aggregate.data.frame(dplyr::select(dd,-date),
                                          by=list(date=agg_period),
                                          FUN=FUN)
    } else {
        if (!is.list(FUN)) {
            stop("FUN should be either a single function or a list of the form list(FUN1=c('var1','var2'), FUN2=c('var3', 'var4'))")
        }
        dd_tmp <- list()
        for (i in seq_along(FUN)) {
            cur_FUN <- get(names(FUN)[i])
            for (j in seq_along(FUN[[i]])) {
                pat <- FUN[[i]][[j]]
                if (!fixed_vars) {
                    var <- grep(pattern=pat,names(dd), value=TRUE)
                } else var <- pat
                if (length(var)==0 || any(!var %in% names(dd))) {
                    warning("no variables matching ",sQuote(pat))
                } else {
                    dd_tmp <- c(dd_tmp,
                                setNames(list(stats::aggregate.data.frame(dd[var],
                                                           by=list(date=agg_period),
                                                           FUN=cur_FUN)[,-1]),var))
                }
            }
            ## FIXME: fix order of columns?
        } ## loop over agg_fun elements
        dd <- do.call(data.frame,c(list(date=unique(na.omit(agg_period))),dd_tmp))
    } ## agg_fun is a list
    dd$date <- as.Date(dd$date)
    dd <- put_attr(dd,aa)
    return(dd)
}

##' summary method for \code{pansim} parameter objects
##'
##' FIXME: add descriptions of what the components are
##' FIXME: add CFR, IFR, etc.
##'
##' @param object a params vector
##' @param ... unused
##' @export
summary.params_pansim <- function(object, ...) {
    check_dots(...)
    ## FIXME: include kappa once we know what works
    ## (analytical vs kernel vs ...)
    if (!"c_prop" %in% names(object)) {
        CFR_gen <- NA
    } else {
        CFR_gen <- with(as.list(object), ((1-alpha)*(1-mu)*(1-phi1)*phi2)/c_prop)
    }
    res <- c(r0=get_r(object),R0=get_R0(object),Gbar=get_Gbar(object),
             CFR_gen=CFR_gen)
    ## FIXME: add IFR, CFR_hosp, CFR_ICU ?
    res["dbl_time"] <- log(2)/res["r0"]
    return(res)
}

##' @export
summary.pansim <- function(object, ...) {
    check_dots(...)
    ## global variables
    ICU <- H <- NULL
    ## FIXME: get ventilators by multiplying ICU by 0.86?
    ## FIXME: prettier?
    p <- attr(object,"params") ## extract params before condensing
    ## test for previous condensation (FIXME)
    if (!"I" %in% names(object)) {
        object <- condense(object)
    }
    unpack(object)
    res <- data.frame(peak_ICU_date=date[which.max(ICU)],
             peak_ICU_val=round(max(ICU)),
             peak_H_date=date[which.max(H)],
             peak_H_val=round(max(H)))
    ## FIXME: report time-varying R0
    if (!is.null(p)) {
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
    psi3 = "1 / mean days post-ICU until discharge",
    ## summary parameters:
    r0 = "initial epidemic growth rate",
    R0 = "basic reproduction number",
    Gbar = "mean generation interval",
    dbl_time = "doubling time",
    ## size arg of rnbinom():
    nb_disp = "dispersion of negative binomial distribution",
    nb_disp.H = "dispersion of negative binomial distribution (hospitalizations)",
    nb_disp.report = "dispersion of negative binomial distribution (case reports)",
    nb_disp.death = "dispersion of negative binomial distribution (deaths)",
    nb_disp.postest = "dispersion of negative binomial distribution (positive tests)",
    zeta = "exponent of phenomenological response to susceptible depletion",
    c_prop = "fraction of infections reported",
    CFR_gen = "case fatality proportion: (1-alpha)*(1-mu)*(1-phi1)*phi2/c_prop",
    W_asymp = "relative testing intensity in asymptomatic compartments"
)

##' Describe parameters
##'
##' Create a data frame with symbols, values and meanings of parameters
##' @param x a \code{params_pansim} object
##' @param stop_missing_names stop if names are missing descriptions? (warn by default)
##' @export
describe_params <- function(x, stop_missing_names=FALSE) {
    if (!is.null(attr(x,"description"))) {
        x_meanings <- attr(x,"description")[names(x)]
    } else {  ## backup/built-in
        m <- match(names(x),names(param_meanings))
        if (any(is.na(m))) {
            wstr <- paste("parameters without description: ",
                          paste(names(x)[is.na(m)],collapse=","))
            if (stop_missing_names) stop(wstr) else warning(wstr)
        }
        x <- x[!is.na(m)]
        x_meanings <- param_meanings[na.omit(m)]
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
##' params <- read_params("ICU1.csv")
##' print(params)
##' print(params,describe=TRUE)
##' @export
## FIXME: prettier printing, e.g. detect "1/" or "proportion"
print.params_pansim <- function( x, describe=FALSE, ... ) {
    check_dots(...)
    if (!describe) {
        attr(x,"description") <- NULL
        print(unclass(x))
    } else {
        print(describe_params(x))
    }
}

##' update parameters within a list of parameters
##' @param object a \code{params_pansim} object
##' @param ... a list of named elements, or a single named list if \code{.list} is specified
##' @param delete_regex vector of regular expressions to remove
##' @param .list treat the first argument as a named list?
##' @export
##' @examples
##' pp <- list(params = c(a=2,b=1))
##' object <- c(a=1, b=0)
##' class(object) <- "params_pansim"
##' update(object,pp,.list=TRUE)
##' pp2 <- list(params = c(a=2))
##' update(object,pp, .list=TRUE)
##' update(object, a=2)
##' update(object, a=2, b=1)
##' update(object, cc1=2, cc2=3)
##' update(object, delete_regex="cc")
update.params_pansim <- function(object, ..., delete_regex=NULL, .list=FALSE) {
    L <- list(...)
    if (.list) {
        if (length(L)>1) stop(".list specified with >1 args")
        L <- L[[1]]
    }
    nm <- names(L)
    for (i in seq_along(L)) {
        ## named sublist
        if (!is.null(nm2 <- names(L[[i]]))) {
            for (j in seq_along(L[[i]])) {
                object[[nm2[j]]] <- L[[i]][j]
            }
        } else {
            ## atomic
            object[[nm[i]]] <- L[[i]]
        }
    }
    if (!is.null(delete_regex)) {
        nm <- names(object)
        for (d in delete_regex) {
            object <- object[!grepl(d,nm)]
        }
    }
    return(object)
}

## need this because of stat/bbmle coef confusion??
##' @method coef fit_pansim
##' @export
coef.fit_pansim <- function(object,
                            method=c("all","fitted"),
                            ...) {
    method <- match.arg(method)
    check_dots(...)
    f_args <- object$forecast_args
    opt_pars <- invlink_trans(restore(coef(object$mle2),
                                      f_args$opt_pars,
                                      f_args$fixed_pars))
    if (method=="fitted") return(opt_pars)
    params <- f_args$base_params
    if (!is.null(opt_pars$params)) params <- update(f_args$base_params, opt_pars$params)
    return(params)
}

##' @export
summary.fit_pansim <- function(object, ...) {
    check_dots(...)
    f_args <- object$forecast_args
    f_args <- legacy_sim_fun(f_args, update=TRUE)
    f_args <- legacy_time_args(f_args, update=TRUE)
    pp <- coef(object,"fitted")
    extra_pars <- pp[!grepl("^params$|nb_disp",names(pp))]
    ## FIXME: get
    time_tab <- with(f_args,sim_fun(params,
                                    extra_pars=extra_pars,
                                    time_args=time_args,
                                    return_timevar=TRUE))
    if ("zeta" %in% names(pp$param)) {
        S_vec <- with(f_args,sim_fun(params,
                                     extra_pars=extra_pars,
                                     time_args=time_args))
    }
    pp_list <- list(coef(object))
    beta0 <- coef(object)[["beta0"]]
    for (i in seq(nrow(time_tab))) {
        pp_list[[i+1]] <- update(pp_list[[1]],
                                 beta0=beta0*time_tab[i,"Relative_value"])
    }
    names(pp_list) <- format(c(f_args$start_date,time_tab$Date))
    ret <- purrr::map_dfr(pp_list,~as.data.frame(rbind(summary(.))),.id="start_date")
    ##  ADD phenom_het stuff here if necessary
    ## browser()
    return(ret)
}


##' @export
update.fit_pansim <- function(object, ...) {
    cc <- object$call
    L <- list(...)
    for (i in seq_along(L)) {
        cc[[names(L)[i]]] <- L[[i]]
    }
    eval.parent(cc)
}

## find call in attribute
## FIXME: DRY?
##' @export
update.pansim <- function(object, ...) {
    cc <- attr(object,"call")
    L <- list(...)
    for (i in seq_along(L)) {
        cc[[names(L)[i]]] <- L[[i]]
    }
    eval.parent(cc)
}

## FIXME: do better than this?
##' @export
print.fit_pansim <- function(x, ...) {
    cat("call:\n")
    print(x$call)
    xx <- x$mle2
    xx@call <- xx@call.orig <- quote(call_deleted)
    print(xx)
}

##' make forecasts from sim
##' @param object a fitted object
##' @param end_date ending date for sim
##' @param stoch stochasticity
##' @param stoch_start stoch starting date
##' @param keep_vars ...
##' @param ensemble run ensemble?
##' @param new_params parameters to update in base parameters (e.g. adding stochastic parameters)
##' @param ... extra args (passed to forecast_ensemble)
##' @param scale_Sigma inflate/deflate covariance matrix
##' @param Sigma covariance matrix
##' @importFrom bbmle coef
##' @export
##' @examples
##' pp1 <- predict(ont_cal1, keep_vars="Rt")
##' ## example of hacking params
##' ont_cal2 <- ont_cal1
##' ont_cal2$forecast_args$base_params["zeta"] <- 4
##' pp2 <- predict(ont_cal2, keep_vars="Rt")
##' ## if zeta is fitted probably need to hack x$mle2@coef, e.g.
##' ont_cal3 <- ont_cal1
##' ## increase beta0 (from -0.34) rather than
##' ## mess with zeta, since phenom het isn't
##' ## estimated in this fit
##' ont_cal3$mle2@fullcoef["params.log_beta0"] <- 0
##' pp3 <- predict(ont_cal3, keep_vars="Rt")
##' pp <- dplyr::bind_rows(base=pp1,zeta=pp2,beta0=pp3, .id="fit")
##' if (require("ggplot2")) {
##'    ggplot(pp,aes(date,value,colour=fit))+geom_line()
##' }
##' \dontrun{
##' ## non-pos-def vcov ... ???
##' predict(ont_cal_2brks,ensemble=TRUE)
##' }
predict.fit_pansim <- function(object
                             , end_date=NULL
                             , stoch=NULL
                             , stoch_start = NULL
                             , keep_vars=c("H","ICU","death", "hosp",
                                           "incidence","report", "cumRep", "newTests/1000")
                             , ensemble = FALSE
                             , new_params=NULL
                             , Sigma=NULL
                             , scale_Sigma=1
                             , ... ) {

    var <- . <- NULL

    f_args <- object$forecast_args
    if ("break_dates" %in% names(f_args)) {
        warning("using old object; switch to time_args(break_dates=...)")
        f_args$time_args <- list(break_dates=f_args$break_dates)
        f_args$break_dates <- NULL
    }
    new_args <- list(...)
    ## FIXME:: dangerous
    for (n in intersect(names(new_args),names(f_args))) {
        f_args[[n]] <- new_args[[n]]
        new_args[[n]] <- NULL
    }
    if (!is.null(end_date)) {
        f_args$end_date <- end_date
    }
    if (!is.null(stoch)) {
        f_args$stoch <- stoch
        f_args$stoch_start <- stoch_start
    }
    if (!is.null(new_params)) {
        f_args$base_params <- do.call(update.params_pansim,
            c(list(f_args$base_params), new_params))
    }
    calc_Rt <- "Rt" %in% keep_vars
    if (!ensemble) {
        fc <- do.call(forecast_sim,
                      c(nlist(p=coef(object$mle2),
                              calc_Rt), f_args, new_args))
    } else {
        ## ensemble
        argList <- c(nlist(fit=object, forecast_args=f_args, scale_Sigma, calc_Rt), new_args)
        ## FIXME: how hard should we work on Sigma?
        if (is.null(Sigma)) {
            if (!is.null(de <- attr(object,"de"))) {
                Sigma <- de$member$Sigma
            } else {
                Sigma <- bbmle::vcov(object$mle2)
            }
            argList <- c(argList,list(Sigma=Sigma))
        }
        fc <- do.call(forecast_ensemble, argList)
    }
    if (inherits(fc,"array")) return(fc)
    fc <- (fc
        %>% sub_vars(keep_vars)
        %>% get_type()
    )
    attr(fc,"forecast_args") <- object$forecast_args
    class(fc) <- c("predict_pansim", class(fc))
    return(fc)
}

## FIXME: don't hard-code!
## data frame for labeling new tests
newtest_lab <-data.frame(date=as.Date("2020-04-10"),
                         value=10,
                         var="newTests/1000",
                         lab="newTests/1000")

## data frame for labeling ICU capacities
capac_info <- data.frame(value=c(630,1300),
                         vtype="prev",
                         var="ICU",
                         lab=c("current","expanded"))

## FIXME: fix upstream/consistently
##  (trans_state_vars?)

##' plot forecasts from fits
##' @param x a calibrated object (result from \code{\link{calibrate}}) or a prediction (from \code{\link{predict.fit_pansim}})
##' @param data original time series data
##' @param break_dates breakpoints
##' @param dlspace spacing for direct labels (not working)
##' @param limspace extra space (in days) to add to make room for direct labels
##' @param add_tests plot newTests/1000?
##' @param add_ICU_cap include horizontal lines showing ICU capacity?
##' @param mult_var variable in data set indicating multiple forecast types to compare
##' @param directlabels use direct labels?
##' @param log use a log10 scale for the y axis?
##' @param log_lwr lower limit when using log scale
##' @param ... extra arguments (unused)
##' @importFrom ggplot2 scale_y_log10 geom_vline facet_wrap theme element_blank geom_line expand_limits geom_point geom_text aes_string labs geom_hline
##' @examples
##' plot(ont_cal1)
##' ont_trans <- trans_state_vars(ont_all)
##' plot(ont_cal1,data=ont_trans)
##' plot(ont_cal1,data=ont_trans, add_tests=TRUE)
##' plot(ont_cal1,data=ont_trans, predict_args=list(end_date="2020-07-01"))
##' \donttest{
##' ## FIXME: don't try these until we have an example where ensemble works
##' ## pp <- predict(ont_cal_2brks, ensemble=TRUE)
##' ## plot(pp)
##' ## plot(pp, data=ont_trans)
##' }
##' @importFrom stats predict
##' @export
plot.predict_pansim <- function(x,
                    data=NULL,
                    break_dates=NULL,
                    dlspace=1,
                    limspace=10,
                    add_tests=FALSE,
                    add_ICU_cap=FALSE,
                    mult_var=NULL,
                    directlabels=TRUE,
                    log=TRUE,
                    log_lwr=1,
                    ...) {
    check_dots(...)
    lwr <- upr <- lab <- var <- . <- NULL
    f_args <- attr(x,"forecast_args")
    if (is.null(break_dates)) break_dates <- legacy_break_dates(f_args)
    var <- date <- value <- mult_var <- NULL

    p <- (ggplot(x,aes(date,value,colour=var))
        + facet_wrap(~vtype,ncol=1,scales="free_y")
        + labs(y="")
        + theme(legend.position="none",
                ## https://stackoverflow.com/questions/10547487/remove-facet-wrap-labels-completely
                strip.background = element_blank(),
                strip.text = element_blank())
    )
    if (log) {
        p <- p + scale_y_log10(limits=c(log_lwr,NA),oob=scales::squish)
    }
    if (!is.null(break_dates)) {
        p <- p + geom_vline(xintercept=as.Date(break_dates),lty=2)
    }
    if (is.null(mult_var)) {
        p <- p + geom_line()
    } else {
        p <- p + geom_line(aes_string(lty=mult_var))
    }
    if (all(c("lwr","upr") %in% names(x))) {
        p <- (p
            + geom_ribbon(aes(ymin=lwr,ymax=upr,fill=var),
                          colour=NA, alpha=0.2)
        )
    }
    if (!is.null(data)) {
        if (add_tests) data <- scale_newtests(data)
        data <- get_type(sub_vars(data))
        p <- (p + geom_point(data=data)
            + geom_line(data=data,alpha=0.2)
        )
        if (add_tests) {
            p <- p + geom_text(data=newtest_lab %>% get_type(),aes(label=lab))
        }
    }
    if (add_ICU_cap) {
        p <- (p
            + geom_hline(data=capac_info,aes(yintercept=value,
                                             colour=var),lty=3)
            + geom_text(data=capac_info,aes(y=value,x=min(data$date),
                                            label=lab),vjust="middle")
        )
    }
    ## trying to fix spacing on the fly: kluge!
    ## dlspace not found
    ## geom_dl() must do weird evaluation env stuff
    if (directlabels) {
        if (!requireNamespace("directlabels")) {
            stop("please install directlabels package")
        }
        p <- (p
            + directlabels::geom_dl(method=list(directlabels::dl.trans(x=x+1),cex=1,'last.bumpup'),
                      aes(label=var))
            + expand_limits(x=max(x$date)+limspace)
        )
    }
    return(p)
}

scale_newtests <- function(x) {
    newTests <- NULL
    xx  <- (x
        %>% dplyr::mutate_at("var",trans_state_vars)
        %>% tidyr::pivot_wider(names_from="var",values_from="value",id_cols="date")
        %>% dplyr::mutate(`newTests/1000`=newTests/1000)
        %>% dplyr::select(-newTests)
        %>% tidyr::pivot_longer(names_to="var",-date)
    )
    return(xx)
}

## this hackery allows us to make plot.fit_pansim a thin layer around
## plot.predict_pansim, with the same arguments (except predict_args,
## which gets used to call predict() and then removed from the argument
## list before calling the predict method for the forecast

##' @rdname predict.fit_pansim
##' @inheritParams plot.predict_pansim
##' @param predict_args additional arguments to pass to predict
##' @export
plot.fit_pansim <- function(x,predict_args=NULL) {
    mc <- match.call()
    forecast <- do.call(predict,c(list(x),predict_args))
    mc[[1]] <- quote(plot)
    mc$x <- forecast
    mc$predict_args <- NULL
    eval.parent(mc)
}

fm0 <- formals(plot.fit_pansim)
fm1 <- formals(plot.predict_pansim)
formals(plot.fit_pansim) <- c(fm0, fm1[-1])

#' @export
vcov.fit_pansim <- function(object, ...) {
    v <- try(solve(object$hessian))
    return(v)
}

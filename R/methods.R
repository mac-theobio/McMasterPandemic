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
                                   c_delay_cv,
                                   max_len=10)
                 )
    ## FIXME: don't hard-code max len ...
    ret <- as.numeric(stats::filter(i,kern,sides=1))
    return(ret)
}

## compute incidence and reports (as convolution of incidence)
calc_reports <- function(x,params) {
    incidence <- x$foi*x$S
    ret <- data.frame(incidence, report=calc_conv(incidence,params))
    return(ret)
}

## FIXME: allow faceting automatically? (each var alone or by groups?)
## don't compare prevalences and incidences?
##' plot method for simulations
##' @param x fitted \code{pansim} object
##' @param drop_states states to \emph{exclude} from plot
##' @param keep_states states to \emph{include} in plot (overrides \code{drop_states})
##' @param condense condense states (e.g. all ICU states -> "ICU") before plotting?  See \code{\link{condense.pansim}}
##' @param log plot y-axis on log scale?
##' @param show_times indicate times when parameters changed?
##' @param ... additional arguments to \code{\link{condense.pansim}}
##' @importFrom ggplot2 ggplot geom_line aes geom_vline scale_y_log10 geom_ribbon
##' @importFrom dplyr one_of
##' @export
plot.pansim <- function(x, drop_states=c("t","S","R","E","I","incidence"),
                        keep_states=NULL, condense=TRUE,
                        log=FALSE, show_times=TRUE, ...) {
    ## global variables
    var <- value <- NULL
    ## attributes get lost somewhere below ...
    ptv <- attr(x,"params_timevar")
    if (condense) {
        x <- condense(x)
    }
    if (!is.null(keep_states)) {
        drop_states <- setdiff(names(x), c(keep_states,"date"))
    }
    ## don't try to drop columns that aren't there
    ## FIXME: use condense.pansim method?
    drop_states <- intersect(drop_states,names(x))
    ## FIXME: don't pivot if already pivoted
    xL <- (pivot(x)
        %>% dplyr::filter(!var %in% drop_states)
        %>% dplyr::mutate(var=forcats::fct_inorder(factor(var)))
        ## FIXME: order factor in pivot?
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

##' Condense columns (infected, ICU, hospitalized) in a pansim output
##' @param object a pansim object
##' @param add_reports add incidence and case reports?
##' @param diff_deaths compute first differences of death series to get daily deaths?
##' @param keep_all keep unaggregated variables in data frame as well?
##' @param ... additional args
##' @export
condense.pansim <-  function(object, add_reports=TRUE, diff_deaths=TRUE, keep_all=FALSE, ...) {
    check_dots(...)
    aa <- get_attr(object)
    ## condense columns and add, if present
    add_col <- function(dd,name,regex) {
        vars <- grep(regex, names(object), value=TRUE)
        if (length(vars)>0) {
            dd[[name]] <- rowSums(object[vars])
        }
        return(dd)
    }
    if (keep_all) {
        dd <- object
    } else {
        dd <- object[c("date","S","E")]
    }
    dd <- add_col(dd,"I","^I[^C]")
    dd <- add_col(dd,"H","^H")
    dd <- add_col(dd,"ICU","^ICU")
    dd <- data.frame(dd,R=object[["R"]])
    dd <- add_col(dd,"discharge","discharge")
    if (diff_deaths) {
        dd <- data.frame(dd,d=c(NA,diff(object[["D"]])))
    } else {
        dd <- data.frame(dd,D=object[["D"]])
    }
    if (add_reports) {
        params <- attr(object,"params")
        if (!"c_delay_mean" %in% names(params)) {
            warning("add_reports requested but delay parameters missing")
        } else {
            dd <- data.frame(dd, calc_reports(object,params))
        }
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
##' sdate <- "10-Feb-2020" ## arbitrary!
##' res <- run_sim(params,state,start_date=sdate,end_date="1-Jun-2020")
##' a1 <- aggregate(res, start="12-Feb-2020",period="7 days", FUN=sum)
##' ## column-specific aggregation
##' first <- dplyr::first
##' a2 <- aggregate(condense(res), start="12-Feb-2020",period="7 days",
##'         FUN=list(mean=c("H","ICU","I"),
##'                sum=c("report","d")))
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
    agg_datevec <- seq.Date(anydate(start),max(dd$date)+extend,
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
##' @param object a params vector
##' @param ... unused
##' @export
summary.params_pansim <- function(object, ...) {
    check_dots(...)
    ## FIXME: include kappa once we know what works
    ## (analytical vs kernel vs ...)
    res <- c(r0=get_r(object),R0=get_R0(object),Gbar=get_Gbar(object))
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
    xa <- condense(object)
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
    psi3 = "1 / mean days post-ICU until discharge",
    ## summary parameters:
    r0 = "initial epidemic growth rate",
    R0 = "basic reproduction number",
    Gbar = "mean generation interval",
    dbl_time = "doubling time"
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
    check_dots(...)
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
        if (length(L[[i]])>1) {
            for (j in seq_along(L[[i]])) {
                object[[names(L[[i]])[j]]] <- L[[i]][j]
            }
        } else {
            ## atomic
            object[[nm[i]]] <- L[[i]]
        }
    }
    return(object)
}

##' @export
summary.fit_pansim <- function(object, ...) {
    check_dots(...)    
    f_args <- attr(object,"forecast_args")
    pars <- invlink_trans(restore(object$par,f_args$opt_pars,f_args$fixed_pars))
    pp <- list()
    beta0 <- pars$params[["beta0"]]
    pp[[1]] <- update(f_args$base_params,beta0=beta0)
    for (i in seq_along(f_args$break_dates)) {
        pp[[i+1]] <- update(pp[[1]],
                            beta0=beta0*pars$rel_beta0[i])
    }
    ## r_jac <- suppressWarnings(vapply(pp,get_r,method="analytic",numeric(1)))
    names(pp) <- c(format(f_args$start_date),format(f_args$break_dates))
    ret <- purrr::map_dfr(pp,~as.data.frame(rbind(summary(.))),.id="start_date")
    return(ret)
}


##' @export
update.fit_pansim <- function(object, ...) {
    cc <- attr(object, "call")
    L <- list(...)
    for (i in seq_along(L)) {
        cc[[names(L)[i]]] <- L[[i]]
    }
    eval.parent(cc)
}

## same strategy: find call in attribute
##' @export
update.pansim <- update.fit_pansim


##' make forecasts from sim
##' @param object a fitted object
##' @param end_date ending date for sim
##' @param stoch stochasticity
##' @param keep_vars ...
##' @param ensemble run ensemble?
##' @param ... extra args (passed to forecast_ensemble)
##' @export
##' @examples
##' predict(ont_cal1)
##' predict(ont_cal_2brks,ensemble=TRUE)
predict.fit_pansim <- function(object
                             , end_date=NULL
                             , stoch=NULL
                             , keep_vars=c("H","ICU","d","incidence","report","newTests/1000")
                              , ensemble = FALSE,
                               ... ) {

    var <- . <- NULL
    get_type <- (.
        %>%  dplyr::mutate(vtype=ifelse(var %in% c("incidence","report","d"),
                                 "inc","prev"))
    )
    sub_vars <- (. %>% dplyr::filter(var %in% keep_vars)
    )
    f_args <- attr(object, "forecast_args")
    if (!is.null(end_date)) {
        f_args$end_date <- end_date
    }
    if (!is.null(stoch)) {
        f_args$stoch <- stoch
    }
    if (!ensemble) {
        fc <- (do.call(forecast_sim,
                       c(list(p=object$par), f_args, list(...))))
    } else {
        argList <- c(list(fit=object, forecast_args=f_args), list(...))
        fc <- do.call(forecast_ensemble, argList)
    }
    fc <- (fc
        %>% sub_vars()
        %>% get_type()
    )
    return(fc)
}

## FIXME: don't hard-code!
## data frame for labeling new tests
newtest_lab <-data.frame(date=as.Date("2020-04-10"),
                         value=10,
                         vtype="prev",
                         var="newTests/1000",
                         lab="newTests/1000")
## data frame for labeling ICU capacities
capac_info <- data.frame(value=c(630,1300),
                         vtype="prev",
                         var="ICU",
                         lab=c("current","expanded"))

##' plot forecasts from fits
##' @param x a calibrated object (result from \code{\link{calibrate}})
##' @param data original time series data
##' @param break_dates breakpoints
##' @param dlspace spacing for direct labels (not working)
##' @param limspace extra space (in days) to add to make room for direct labels
##' @param add_tests plot newTests/1000?
##' @param predict_args arguments to pass to predict method
##' @param add_ICU_cap include horizontal lines showing ICU capacity?
##' @param mult_var variable in data set indicating multiple forecast types to compare
##' @param directlabels use direct labels?
##' @param ... extra arguments (unused)
##' @importFrom ggplot2 scale_y_log10 geom_vline facet_wrap theme element_blank geom_line expand_limits geom_point geom_text aes_string labs geom_hline
##' @importFrom directlabels geom_dl dl.trans
##' @examples
##' plot(ont_cal1)
##' ont_trans <- trans_state_vars(ont_all)
##' plot(ont_cal1,data=ont_trans)
##' plot(ont_cal1,data=ont_trans, add_tests=TRUE)
##' plot(ont_cal1,data=ont_trans, predict_args=list(end_date="2020-07-01"))
##' \donttest{
##' plot(ont_cal_2brks,predict_args=list(ensemble=TRUE))
##' }
##' @importFrom stats predict
##' @export
plot.fit_pansim <- function(x,
                    data=NULL,
                    break_dates=NULL,
                    dlspace=1,
                    limspace=10,
                    add_tests=FALSE,
                    add_ICU_cap=FALSE,
                    mult_var=NULL,
                    predict_args = NULL,
                    directlabels=TRUE,
                    ...) {
    check_dots(...)
    lwr <- upr <- lab <- var <- . <- NULL
    f_args <- attr(x,"forecast_args")
    if (is.null(break_dates)) break_dates <- f_args$break_dates
    forecast <- do.call(predict,c(list(x),predict_args))
    var <- date <- value <- mult_var <- NULL
    p <- (ggplot(forecast,aes(date,value,colour=var))
        + scale_y_log10(limits=c(1,NA),oob=scales::squish)
        + facet_wrap(~vtype,ncol=1,scales="free_y")
        + labs(y="")
        + theme(legend.position="none",
                ## https://stackoverflow.com/questions/10547487/remove-facet-wrap-labels-completely
                strip.background = element_blank(),
                strip.text = element_blank())
    )
    if (!is.null(break_dates)) {
        p <- p + geom_vline(xintercept=break_dates,lty=2)
    }
    if (is.null(mult_var)) {
        p <- p + geom_line()
    } else {
        p <- p + geom_line(aes_string(lty=mult_var))
    }
    if (all(c("lwr","upr") %in% names(forecast))) {
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
            p <- p + geom_text(data=newtest_lab,aes(label=lab))
        }
    }
    if (add_ICU_cap) {
        p <- (p
            + geom_hline(data=capac_info,aes(yintercept=value,
                                             colour=var),lty=3)
            + geom_text(data=capac_info,aes(y=value,x=min(data$date),
                                            label=lab),vjust=-1)
        )
    }
    ## trying to fix spacing on the fly: kluge!
    ## dlspace not found
    ## geom_dl() must do weird evaluation env stuff
    if (directlabels) {
        p <- (p
            + geom_dl(method=list(dl.trans(x=x+1),cex=1,'last.bumpup'),
                      aes(label=var))
            + expand_limits(x=max(forecast$date)+limspace)
        )
    }
    return(p)
}

## FIXME: less hard-coding
keep_vars <- c("H","ICU","d","incidence","report","newTests/1000")

get_type <- . %>%  dplyr::mutate(vtype=ifelse(var %in% c("incidence","report","d"),
                                       "inc","prev"))
sub_vars <- . %>% dplyr::filter(var %in% keep_vars)


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

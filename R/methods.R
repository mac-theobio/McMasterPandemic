##' @export
print.pansim <- function(x,all=FALSE,...) {
    ## FIXME: is this the best way?
    ## use tibbles or not?
    class(x) <- "data.frame"
    if (all) print(x)
    attr(x,"params") <- NULL
    print(x)
}

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
        %>% dplyr::select(-one_of(drop_states))
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

##' Collapse columns (infected, ICU, hospitalized) in a pansim output
##' @param x a pansim object
##' @param pivot return long-format tibble instead of wide data frame?
##' @param keep_vars variables to retain (in addition to date) if pivoting
##' @param ... unused, for generic consistency
##' @importFrom stats aggregate
##' @importFrom dplyr %>% as_tibble
##' @importFrom tidyr pivot_longer
##' @export
aggregate.pansim <- function(x,pivot=FALSE,keep_vars=c("H","ICU","D"), ...) {
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
    if (!pivot) return(dd)
    ## more convenient for regressions etc.
    dd <- (dd
        %>% as_tibble()
        %>% dplyr::select(c("date",keep_vars))
        %>% pivot_longer(names_to="var",-date)
    )
    return(dd)
}

##' @export
summary.panparams <- function(object, ...) {
    c(r=get_r(object),R0=get_R0(object),get_GI_moments(object))
}

##' @export
summary.pansim <- function(object, ...) {
    ## FIXME: get ventilators by multiplying ICU by 0.86?
    ## FIXME: prettier?
    xa <- aggregate(object)
    attach(xa); on.exit(detach(xa))
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

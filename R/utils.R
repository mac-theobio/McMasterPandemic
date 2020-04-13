
##' attempt to convert x to a date unless it already is one
##' @param x a date, or a character string in some day-month-year format that \code{lubridate::dmy} can handle
##' @keywords internal
##' @export
ldmy <- function(x) if (inherits(x,"Date")) x else lubridate::dmy(x)

##' self-naming list (copied from lme4:::namedList)
##' @param ... a list of objects
##' @export
nlist <- function (...) {
    L <- list(...)
    snm <- vapply(substitute(list(...)), deparse, character(1))[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

## select every ndt'th row of a data frame
thin <- function(x,ndt=1) {
    if(ndt==1) return(x)
    x  <- x[seq(nrow(x)) %% ndt == 1,]
    return(x)
}

##' unpack a list into the current (calling) environment
##' this is a replacement for \code{with()} (which is hard to debug) and \code{attach()} (which is frowned upon by CRAN/triggers package-check warnings)
##' (for more traditional Python tuple unpacking see the \code{zeallot} package)
##' @param x a named vector or a list
##' @export
unpack <- function(x) {
    invisible(list2env(as.list(x),envir=parent.frame()))
}

## generate " a <- b <- c <- NULL" statements
if (FALSE) {
    p <- read_params(system.file("params","ICU1.csv",
                                 package="McMasterPandemic"))
    s <- make_state(params=p)
    m <- make_ratemat(s,p)
    print_globals(s,p)
}
print_globals <- function(..., chunksize=8) {
    all_names <- unlist(lapply(list(...),
              function(z) if (is.matrix(z)) colnames(z) else names(z)))
    ss <- split(all_names,seq(all_names) %/% chunksize)
    ss2 <- paste(sapply(ss,paste,collapse=" <- ")," <- NULL")
    cat(ss2,sep="\n")
}

## utility for Stanford-to-McMaster notation change
rename_params <- function(p) {
    n <- names(p)
    n <- gsub("^gamma","sigma",n)
    n <- gsub("^lambda","gamma",n)
    names(p) <- n
    return(p)
}



## dictionary; internal name, graph label
label_dict <- read.csv(stringsAsFactors=FALSE,
text="
Symbol,Label,Regex
S,Susceptible,Sus.*
E,Exposed,Exp.*
I,Total infectious,Inf.*
Ia,Infectious/asymptomatic,NA
Im,Infectious/mild,NA
Is,Infectious/severe,NA
H,Hospital,Hosp.*
ICU,ICU,ICU.*
D,Deaths,[Dd]e.*
R,Recovered,NA
report,Case reports,newConf.*
")
## factor with levels in order of appearance (forcats equiv?)
ff <- function(x) factor(x,levels=x)
label_dict$Symbol <- ff(label_dict$Symbol)
label_dict$Label <-  ff(label_dict$Label)

##' Translate state variables
## FIXME: NA for non-matches??? cleaner way to do this?
##' @importFrom stats na.omit
##' @param x a vector of state variables (factor or character)
##' @export
trans_state_vars <- function(x) {
    if (inherits(x,"data.frame")) stop("trans_state_vars should be applied to a vector, not a data frame")
    x <- as.character(x)
    D <- na.omit(label_dict)
    matches <- lapply(D$Regex,grep,x=unique(x))
    if (any(vapply(matches,length,numeric(1))>1)) {
        stop("multiple matches")
    }
    for (i in seq(nrow(D))) {
        if (length(matches[[i]])>0) {
            ## translate
            x[grepl(D$Regex[i],x)] <- as.character(D$Symbol[i])
        }
    }
    return(x)
}        

##' read simulation parameters
##' @param fn file name (CSV file containing at least value and symbol columns);  file should either be findable from current working directory or available in the built-in \code{params} directory
##' @param value_col name of column containing values
##' @param symbol_col name of column containing symbols
##' @param desc_col name of (optional) column containing descriptions
##' @importFrom stats setNames
##' @importFrom utils read.csv write.table
##' @export
read_params <- function(fn,value_col="Value",symbol_col="Symbol",
                        desc_col="Parameter") {
    if (!file.exists(fn)) {
        fn <- system.file("params",fn,package="McMasterPandemic")
    }
    x <- read.csv(fn,
                  colClasses="character",
                  stringsAsFactors=FALSE,
                  comment.char="#",
                  na.strings="variable")
    ## evaluate to allow expressions like "1/7" -> numeric
    x[[value_col]] <- vapply(x[[value_col]], function(z) eval(parse(text=z)), numeric(1))
    res <- setNames(x[[value_col]],x[[symbol_col]])
    class(res) <- "params_pansim"
    if (desc_col %in% names(x)) {
        attr(res,"description") <- setNames(x[[desc_col]],x[[symbol_col]])
    }
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

##' apply inverse-link functions to parameter vector or list based on names
##' @param p a named vector or list of parameter values with names of the form <linkfun>_name
##' @examples
##' invlink_trans(c(log_p1=0,logit_p2=0))
##' invlink_trans(list(log_p1=c(0,0),logit_p2=c(0,0,0)))
##' invlink_trans(list(p1=c(log_a=0,log_b=0),p2=4))
##' @export
invlink_trans <- function(p) {
    r <- vector("list",length(p))
    for (i in seq_along(p)) {
        ## recurse if necessary
        if (length(p[[i]])>1 && !is.null(names(p[[i]]))) {
            r[[i]] <- invlink_trans(p[[i]])
        } else {
            nm <- names(p)[i]
            if (!grepl("_",nm)) {
                ## identity; should be able to do this with a better regex?
                r[[i]] <- p[[i]]
            } else {
                invlink <- gsub("^([^_]+).*","\\1",names(p)[i])
                ## cat(invlink,"\n")
                r[[i]] <- switch(invlink,
                                 log=exp(p[[i]]),
                                 log10=10^(p[[i]]),
                                 logit=plogis(p[[i]]),
                                 stop("unknown link"))
                ## FIXME: add cloglog? user-specified links?
            } ## contains link
        }  ## atomic
    } ## loop over p
    names(r) <- gsub("^[^_]+_","",names(p))
    if (is.numeric(p)) r <- unlist(r)
    return(r)
}

##' convert a parameter vector back to a structured list
##' Like \code{\link{relist}}, but adds the capability to fill in
##' parameters from a "fixed-parameter" vector
##' @param flesh a vector to be restored
##' @param skeleton a list, the structure of which determines the structure of the result
##' @param fixed a list which determines extra components to fill in
##' @note Depends at present on the names of the unlisted object; may be fragile
##' @examples
##' opt_pars <- list(log_E0=4, log_beta0=-1, log_rel_beta0=c(-1,-1), log_nb_disp=0)
##' restore(unlist(opt_pars),opt_pars)
##' invlink_trans(restore(unlist(opt_pars),opt_pars))
##' @export
restore <- function(flesh, skeleton, fixed=NULL) {
    flesh <- unlist(flesh)  ## just in case ...
    if (is.null(fixed)) return(utils::relist(flesh,skeleton))
    ## rely on name-matching for now ... fragile?
    full_flesh <- unlist(skeleton)
    full_flesh[names(flesh)] <- flesh
    full_flesh[names(fixed)] <- unlist(fixed)
    utils::relist(full_flesh,skeleton)
}

## save attributes
get_attr <- function(x, exclude=c("names", "row.names", "class")) {
    aa <- attributes(x)
    aa <- aa[setdiff(names(aa), exclude)]
    return(aa)
}

## restore attributes and set class to (pansim, ...)
put_attr <- function(x, a) {
    for (i in seq_along(a)) {
        attr(x,names(a)[i]) <- a[i]
    }
    class(x) <- c("pansim", class(x))
    return(x)
}

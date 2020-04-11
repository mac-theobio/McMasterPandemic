
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

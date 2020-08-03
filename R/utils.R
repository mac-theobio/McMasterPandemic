## FIXME: less hard-coding!!!
keep_vars <- c("H","ICU","death", "incidence","report","newTests/1000", "hosp")
get_type <- . %>%  dplyr::mutate(vtype=ifelse(var %in% c("incidence","report","death","newTests/1000", "hosp"),  "inc","prev"))
## un-tidyverse this because I don't want to learn how to do NSE right now
sub_vars <- function(x,kv=keep_vars) { if (identical(kv,"all")) x else x[x$var %in% kv,]  }

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
## general convention: capital letters for prevalences or cumulative values
## lowercase for incidences
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
D,Deaths,^(D|dec)
death,New deaths,newDeath.*
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
    if (inherits(x,"data.frame")) {
        return(dplyr::mutate_at(x,"var",trans_state_vars))
    }
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
##' @param p a named vector or list of parameter values with names of the form \code{linkfun_name}
##' @param unknown_link action to take when unknown link is specified (ignored)
##' @param links allowed link prefixes
##' @examples
##' invlink_trans(c(log_p1=0,logit_p2=0))
##' invlink_trans(list(log_p1=c(0,0),logit_p2=c(0,0,0)))
##' invlink_trans(list(p1=c(log_a=0,log_b=0),p2=4))
##' tst <- list(params = c(log_beta0 = 0.693147180559945), log_nb_disp = 4.60517018598809)
##' invlink_trans(tst)
##' tst2 <- list(params=c(log_E0=1,log_beta0=1),log_nb_disp=c(H=0,report=0,death=0))
##' invlink_trans(tst2)
##' invlink_trans(list(time_beta=1,log_time_beta=0))  ## ignore non-link prefixes
##' invlink_trans(list(time_beta=c(2,3),log_time_beta=c(a=0,b=0)))  ## 
##' @export
invlink_trans <- function(p, unknown_link="ignore", links=c("log","log10","logit")) {
    r <- vector("list",length(p))
    for (i in seq_along(p)) {
        nm <- names(p)[i]
        invlink <- gsub("^([^_]+)_.*","\\1",nm)
        has_link <- invlink %in% links
        ## recurse if necessary
        ## if ((length(p[[i]])>1 || is.list(p[[i]])) && !is.null(names(p[[i]]))) {
        if (!is.null(names(p[[i]]))) {
            pp <- p[[i]]
            ## attach link name to sub-parameters
            if (has_link) pp <- setNames(pp,paste(invlink,names(p[[i]]),sep="_"))
            r[[i]] <- invlink_trans(pp)
            if (has_link) names(r)[[i]] <- gsub("^[^_]+_","",nm) else names(r)[[i]] <- nm ## ?? not sure
        } else {
            if (!has_link) {
                ## identity; should be able to do this with a better regex?
                r[[i]] <- p[[i]]
                names(r)[i] <- nm
            } else {
                ## cat(invlink,"\n")
                r[[i]] <- switch(invlink,
                                 log=exp(p[[i]]),
                                 log10=10^(p[[i]]),
                                 logit=plogis(p[[i]]),
                                 { if (unknown_link=="stop") stop("unknown link: %s",invlink)
                                     if (unknown_link=="warn") warning("unknown link: %s",invlink)
                                     p[[i]]
                                 })
                names(r)[i] <- gsub("^[^_]+_","",nm)
                ## FIXME: add cloglog? user-specified links
            } ## contains link
        }  ## atomic
    } ## loop over p
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
## FIXME: fixed param might be obsolete (and we could go back to just using relist()
## now that we're using mle2 with its own fixed capabilities?
restore <- function(flesh, skeleton, fixed=NULL) {
    if (length(nf <- names(flesh))>0) {
        ## FIXME: add test here for matching names?
    }
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

## action: message, warning, stop
check_dots <- function(..., action="stop") {
    L <- list(...)
    if (length(L)>0) {
        FUN <- get(action)
        FUN("unknown arguments: ",
            paste(names(L), collapse=","))
    }
    return(NULL)
}

#' Convert expression string to TeX format
#'
#' This has very limited capabilities and is intended
#' mainly to convert parameter names in mathematical expression strings
#' to associated TeX symbols so they look good when \code{\link[tikzDevice:tikzDevice-package]{tikz}}
#' is in use.
#'
#' @param x character string to be interpreted
#' @param dollars if \code{TRUE} then surround output string in \code{$} signs
#' @param force if \code{TRUE} then convert to TeX even if \code{\link{dev_is_tikz}} returns \code{FALSE}
#' @examples
#' texify("R0 = beta/gamma")
#' texify("R0 = beta/gamma", dollars=FALSE, force=TRUE)
#' 
#' @seealso \code{\link[Hmisc:latex]{latexTranslate}}
#' @export
texify <- function( x, dollars=TRUE, force=dev_is_tikz() ) {
  if (!is.character(x)) return(x)
  if (!force) {
    ##x <- gsub("fun","",x,fixed=TRUE)
    return(x)
  }
  x <- gsub("E0","E(0)",x,fixed=TRUE)
  x <- gsub("R0","{\\mathcal R}_0",x,fixed=TRUE)
  x <- gsub("r0","r_0",x,fixed=TRUE)
  x <- gsub("Gbar","\\bar{G}",x,fixed=TRUE)
  x <- gsub("dbl_time","T_2",x,fixed=TRUE)
  x <- gsub("beta0","beta_0",x,fixed=TRUE)
  x <- gsub("phi1","phi_1",x,fixed=TRUE)
  x <- gsub("nb_disp","n",x,fixed=TRUE)
  x <- gsub("n.H","n_{\\rm H}",x,fixed=TRUE)
  x <- gsub("n.report","n_{\\rm report}",x,fixed=TRUE)
  x <- gsub("n.death","n_{\\rm death}",x,fixed=TRUE)
  x <- gsub("Ca","C_{\\rm a}",x,fixed=TRUE)
  x <- gsub("Cp","C_{\\rm p}",x,fixed=TRUE)
  x <- gsub("Cs","C_{\\rm s}",x,fixed=TRUE)
  x <- gsub("Cm","C_{\\rm m}",x,fixed=TRUE)
  x <- gsub("gamma_a","gamma_{\\rm a}",x,fixed=TRUE)
  x <- gsub("gamma_p","gamma_{\\rm p}",x,fixed=TRUE)
  x <- gsub("gamma_s","gamma_{\\rm s}",x,fixed=TRUE)
  x <- gsub("gamma_m","gamma_{\\rm m}",x,fixed=TRUE)
  ##FIX: the following would be better using a table of
  ##     Greek letters and/or a table of TeX symbols
  x <- gsub("alpha","\\alpha",x,fixed=TRUE)
  x <- gsub("beta","\\beta",x,fixed=TRUE)
  x <- gsub("gamma","\\gamma",x,fixed=TRUE)
  x <- gsub("delta","\\delta",x,fixed=TRUE)
  x <- gsub("epsilon","\\varepsilon",x,fixed=TRUE)
  x <- gsub("rho","\\rho",x,fixed=TRUE)
  x <- gsub("sigma","\\sigma",x,fixed=TRUE)
  x <- gsub("nu","\\nu",x,fixed=TRUE)
  x <- gsub("mu","\\mu",x,fixed=TRUE)
  x <- gsub("zeta","\\zeta",x,fixed=TRUE)
  x <- gsub("phi","\\phi",x,fixed=TRUE)
  x <- gsub("c_prop","c_{\\rm prop}",x,fixed=TRUE)
  if (dollars) x <- paste0("$",x,"$")
  return(x)
}

#' Is \code{\link[tikzDevice:tikzDevice-package]{tikzDevice}} currently in use?
#'
#' Convenient for dealing with text in graphics,
#' which can be rendered much more professionally with \code{\link[tikzDevice:tikzDevice-package]{tikz}}.
#' @return logical
#' @export
dev_is_tikz <- function() {
  return(names(grDevices::dev.cur()) == "tikz output")
}

## especially for 3.6/4.0 compatibility
dfs <- function(...) data.frame(..., stringsAsFactors=FALSE)

legacy <- function(x, nm, other) {
    if (nm %in% names(x)) x[[nm]] else other
}
## support switch in break_dates specification
legacy_break_dates <- function(x, update=FALSE) {
    r <- legacy(x, "break_dates", x$time_args$break_dates)
    if (!update) return(r)
    x$time_args$break_date <- r
    return(x)
}

## support switch in break_dates specification
legacy_sim_fun <- function(x, update=FALSE) {
    r <- legacy(x, "sim_fun", run_sim_break)
    if (!update) return(r)
    x$sim_fun <- r
    return(x)
}

legacy_time_args <- function(x, update=FALSE) {
    r <- legacy(x, "time_args", list(break_dates=legacy_break_dates(x)))
    if (!update) return(r)
    x$time_args <- r
    return(x)
}


##' generate reasonable default opt_pars
##' @param params baseline parameters
##' @param vars which variables are being fitted to?
##' @examples
##' params <- read_params("ICU1.csv")
##' get_opt_pars(params)
##' get_opt_pars(params,vars="hosp")
##' get_opt_pars(params,vars="report")
##' @importFrom stats qlogis
##' @export
get_opt_pars <- function(params,vars=c("hosp","death","report")) {
    opt_pars <- list(params =  c(log_E0=log(params[["E0"]])
                               , log_beta0=log(params[["beta0"]])))
    if ("death" %in% vars && any(c("hosp","H") %in% vars)) {
        ## fraction of hosp to acute (non-ICU) (\propto death)
        opt_pars$params <- c(opt_pars$params,
                             logit_phi1=qlogis(params[["phi1"]]))
    }
    ## fraction of mild (non-hosp) cases (\propto hosp)
    if (any(c("hosp","H") %in% vars)) {
        opt_pars$params <- c(opt_pars$params, log_mu=log(params[["mu"]]))
    }
    if (identical(sort(vars), c("death","report"))) {
        opt_pars$params <- c(opt_pars$params, logit_c_prop=qlogis(params[["c_prop"]]))
    }
    ## per-parameter dispersion
    opt_pars$log_nb_disp <- setNames(rep(0,length(vars)),vars)
    return(opt_pars)
}

get_DE_lims <- function(opt_pars,default=c(lwr=-1,upr=1),
                        special=list(lwr=c("params.log_E0"=1,zeta=-2,
                                           "time_beta"=-2),
                                     upr=c("rel_beta0"=4,
                                           "nb_disp|E0"=5,
                                           "zeta"=5,
                                           "time_beta"=2))) {
    opt_inputs <- unlist(opt_pars)
    lwr <- opt_inputs  ## get all names
    lwr[] <- default[["lwr"]]
    for (i in seq_along(S <- special[["lwr"]])) {
        lwr[grepl(names(S)[i],names(lwr))] <- S[[i]]
    }
    upr <- opt_inputs  ## get all names
    upr[] <- default[["upr"]]
    for (i in seq_along(S <- special[["upr"]])) {
        upr[grepl(names(S)[i],names(upr))] <- S[[i]]
    }
    return(list(lwr=lwr,upr=upr))
}

##' recursively log-ify expressions
##' @export
##' @param x an expression
##' @examples
##' add_d_log(~dnorm(a,b,c))
##' add_d_log(~sum(dnorm(a,b,c)))
##' @keywords internal
add_d_log <- function(x) {
    if (length(x)==1) return(x)
    nm <- deparse(x[[1]])
    if (substr(nm,1,1)=="d" && nm %in% ls("package:stats")) {
        x$log <- TRUE
    }
    for (i in 2:length(x)) {
        x[[i]] <- add_d_log(x[[i]])
    }
    return(x)
}

##' extract vector of effective Rt
## FIXME: now redundant with code inside forecast_sim
## if this gets erased, save the @importFrom stuff!
##' @param x a fitted object
##' @importFrom dplyr rename mutate_at filter transmute mutate arrange full_join as_tibble case_when
##' @export
get_Rt <- function(x) {
    Symbol <- value <- rel_beta0 <- S <- zeta <- NULL
    ## process args etc.
    f_args <- x$forecast_args
    pp <- coef(x,"fitted")
    extra_pars <- pp[!grepl("^params$|nb_disp",names(pp))]
    ## calc R0 base, rel beta, S
    R0_base <- summary(coef(x))[["R0"]]
    bb <- with(f_args,sim_fun(params=coef(x),
                              extra_pars=extra_pars,
                              time_args=time_args,sim_args=sim_args,
                              return_timevar=TRUE))
    S_pred <- (predict(x,keep_vars="S")
        %>% select(date,value)
        %>% rename(S="value")
    )
    ## combine rel_beta and S information
    if (!is.null(bb)) {
        bb <- (bb
            %>% as_tibble()
            %>% dplyr::filter(Symbol=="beta0")
            %>% dplyr::select(-Symbol)
            %>% rename(rel_beta0="Relative_value",date="Date")
        )
        x2 <- full_join(bb,S_pred,by="date") %>% arrange(date)
    }  else {
        x2 <- mutate(S_pred,rel_beta0=1)
    }
    ## fill in values preceding/following estimated beta0 values
    x2_nona <- na.omit(x2)
    x3 <- (x2 
        %>% mutate_at("rel_beta0",
                      ~ case_when(
    is.na(.) & date<x2_nona$date[1] ~ x2_nona$rel_beta0[1],
    is.na(.) & date>last(x2_nona$date) ~ last(x2_nona$rel_beta0),
    TRUE ~ .
)))
    if (has_zeta(coef(x))) {
        N <- coef(x)[["N"]]
        ## FIXME: use hetS!
        x4 <- (x3
            %>% transmute(date=date,
                          R0=R0_base*rel_beta0*(S/N)^(1+zeta))
        )
    }
    return(x4)
}
    

## replace values before first non-NA value with first non-NA value; na.locf everything else
fill_edge_values <- function(x) {
    first_val <- which(!is.na(x))[1]
    if (first_val>1) {
        x[seq(first_val-1)] <- x[first_val]
    }
    x <- zoo::na.locf(x)
    return(x)
}

## save file with better compression
## ff <- list.files(pattern="\\.rda$")
## lapply(ff, recompress)
recompress <- function(fn) {
    L <- load(fn)
    invisible(save(list=L, file=fn, compress="xz"))
}



#' extract and return original fitting data
#' @param x a fitted object
#' @export
getData <- function(x) {
    x$mle2@data$data
}

has_zeta <- function(params) {
    "zeta" %in% names(params) && params[["zeta"]]!=0
}

## test based *either* on state or params
## testing based on params fails when we have make_state -> get_evec -> make_state ...
has_testing <- function(state,params=NULL) {
    if (missing(state)) {
        return("testing_intensity" %in% names(params) && params[["testing_intensity"]]!=0)
    } else {
        return(any(grepl("_t$",names(state))))
    }
}

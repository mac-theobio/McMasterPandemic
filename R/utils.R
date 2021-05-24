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
    if (any(names(x)=="")) stop("unnamed elements in x")
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
  x <- gsub("phi2","phi_2",x,fixed=TRUE)
  x <- gsub("psi1","psi_1",x,fixed=TRUE)
  x <- gsub("psi2","psi_2",x,fixed=TRUE)
  x <- gsub("psi3","psi_3",x,fixed=TRUE)
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
  x <- gsub("W_asymp","W_{\\rm asymp}",x,fixed=TRUE)
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
  x <- gsub("psi","\\psi",x,fixed=TRUE)
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

## set limits for DEoptim hypercube
##  names of 'special' args are REGULAR EXPRESSIONS
##  need to match parameters as they appear in unlist(opt_pars):
##  values in parameter vector: params.[LINKFUN]_NAME
##  dispersion parameters: log_nb_disp.VAR
##  time parameters: time_beta.???
##  rel_beta0 is from breakpoints models
##  FIXME: can this go in an input file instead?
get_DE_lims <- function(opt_pars,default=c(lwr=-1,upr=1),
                        special=list(lwr=c(log_E0=1,
                                           zeta=-2,
                                           time_beta=-2,
                                           log_testing_intensity=-5,
                                           mu=-1),
                                     upr=c(rel_beta0=4,
                                           nb_disp=5,
                                           E0=5,
                                           zeta=5,
                                           time_beta=2,
                                           log_testing_intensity=-2,
                                           mu=3)))
{
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

## extract vector of effective Rt
## now redundant with code inside forecast_sim
## fix this up as a standalone function, then replace stuff inside forecast_sim
## @param x a fitted object
##'@importFrom dplyr rename mutate_at filter transmute mutate arrange full_join as_tibble case_when
## @export
## get_Rt <- function(x, new_params=NULL, S=NULL) {
##     f_args$base_params <- do.call(update.params_pansim,
##                                   c(list(f_args$base_params), new_params))

##     R0_base <- summary(params)[["R0"]]
##         ## retrieve time-varying beta
##         bb <- do.call(sim_fun,c(all_sim_args,list(return_timevar=TRUE)))
##         if (!is.null(bb)) {
##             bb <- (bb
##                 %>% as_tibble()
##                 %>% dplyr::filter(Symbol=="beta0")
##                 %>% dplyr::select(-Symbol)
##                 %>% rename(rel_beta0="Relative_value",date="Date")
##             )
##             vars <- c("date","S")
##             if (has_zeta(params)) vars <- c(vars,"hetS")
##             x2 <- (full_join(bb,select(r_agg,vars),
##                              by="date")
##                 %>% arrange(date))
##         }  else {
##             x2 <- r_agg %>% select(date,S, hetS) %>% mutate(rel_beta0=1)
##         }
##         x3 <- (x2
##             %>% mutate_at("rel_beta0", fill_edge_values)
##             %>% transmute(date=date, hetS=hetS, Rt=R0_base*rel_beta0)
##         )
##         if (has_zeta(params)) {
##             x3 <- (x3
##                 %>% mutate_at("Rt", ~.*hetS)
##                 %>% select(-hetS)
##             )
##         }
##         r_agg <- full_join(r_agg, x3, by="date")

## }


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

has_vacc <- function(params) {
    "vacc" %in% names(params)
}

## test based *either* on state or params
## testing based on params fails when we have make_state -> get_evec -> make_state ...
has_testing <- function(state,params=NULL,ratemat=NULL) {
    if (!is.null(params)) {
        return("testing_intensity" %in% names(params) && params[["testing_intensity"]]!=0)
    }
    if (!is.null(ratemat)) {
        return(any(grepl("_t$",rownames(ratemat))))
    }
    return(any(grepl("_t$",names(state))))
}

has_age <- function(x) {
  ## look for presence of the "age_cat" attribute
  if("pansim" %in% class(x)){
    ## for pansim objects, check if its params attribute has age cats
    return("age_cat" %in% names(attributes(attr(x, "params"))))
  }
  ## otherwise, check the object directly for an age cat attribute
  return("age_cat" %in% names(attributes(x)))
}

get_age <- function(x) {
  ## get age categories out of params list
  if (!has_age(x)) stop("these parameters are not age-specific")
  return(attr(x, "age_cat"))
}

has_vax <- function(x) {
  ## look for presence of the "vax_cat" attribute
  if("pansim" %in% class(x)){
    ## for pansim objects, check if its params attribute has age cats
    return("vax_cat" %in% names(attributes(attr(x, "params"))))
  }
  ## otherwise, check the object directly for an age cat attribute
  return("vax_cat" %in% names(attributes(x)))
}

get_vax <- function(x) {
  ## get age categories out of params list
  if (!has_vax(x)) stop("these parameters are not vaxified")
  return(attr(x, "vax_cat"))
}

## round, preserving sum
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
  y <- floor(x)
  indices <- utils::tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}

##' visualize rate (per-capita flow) matrix
##' @param M rate matrix
##' @param method visualization method
##' @param subset list of two regular expressions, the first to subset rate matrix rows (based on rownames) and the second to subset columns (based on colnames)
##' @param aspect aspect ratio ("iso", "fill" are the sensible options)
##' @param block_size numeric vector of number of compartments per block; if NA, try to guess from number of epidemiological compartments
##' @param block_col (numeric vector, of length 1 or length(block_size)
##' @param const_width set flows to constant value of 1?
##' @param colour_palette vector of colours for rate matrix heatmap
##' @param do_symbols plot symbolic values for flows?
##' @param axlabs for flow matrices, show axis tick labels?
##' @param box.size box size for diagram
##' @param block_col each element in block_col controls the color of the corresponding grid overlay added by block_size (if add_blocks==TRUE)
##' @param ... arguments to pass to lower level functions (plotmat::diagram/image/igraph)
##' @importFrom lattice panel.abline
##' @importFrom Matrix Matrix
##' @importFrom graphics image
##' @importFrom diagram plotmat
## See \code{help("Matrix::image-methods")} for more.
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params[["N"]],E0=params[["E0"]], use_eigvec=FALSE)
##' M <- make_ratemat(state, params)
##' show_ratemat(M)
##' ## silly but shows we can do multiple block types in different colours
##' show_ratemat(M, block_size=c(3,5), block_col=c(2,4))
##' @export
show_ratemat <- function(M, method=c("Matrix","diagram","igraph"),
                         subset = NULL,
                         xlab = "to",
                         ylab = "from",
                         sub = "",
                         zlim = c(0,1),
                         aspect="iso",
                         block_size=NULL,
                         block_col=2,
                         axlabs=TRUE,
                         const_width=(method=="igraph"),
                         colour_palette=viridis::magma(n=50, direction = -1),
                         do_symbols=NULL,
                         box.size=0.02,...) {
    method <- match.arg(method)
    ## subset ratemat, if desired
    if(!is.null(subset)){
      M <- M[grepl(subset[1], dimnames(M)$from),
             grepl(subset[2], dimnames(M)$to)]
    }

    p <- NULL
    if (is.null(do_symbols)) {
        do_symbols <- method=="diagram" && !has_testing(ratemat=M)
    }
    if (const_width && !do_symbols) { M[M>0] <- 1 }
    if (method=="Matrix") {
        add_blocks <- !is.null(block_size)
        if (is.null(add_blocks)) {
            add_blocks <- has_testing(state=setNames(numeric(nrow(M)),
                                                     rownames(M)))
        }
        if (axlabs) {
            rlabs <- colnames(M)
            clabs <- rownames(M)
        } else {
            rlabs <- clabs <- rep("",nrow(M))
        }
        p <- Matrix::image(Matrix(M),
                           scales=list(x=list(at=seq(ncol(M)),labels=rlabs,
                                              rot=90),
                                       y=list(at=seq(nrow(M)),labels=clabs)),
                           xlab=xlab,
                           ylab=ylab,
                           sub=sub,
                           colorkey = !const_width,
                           col.regions = colour_palette,
                           aspect=aspect, ...)
         if (add_blocks) {
             if (all(is.na(block_size))) {
               ## FIXME: test more. Works suboptimally for a *single* block, but ???
                 epi_vars <- unique(gsub("_.*$","",colnames(M)))
                 epi_vars <- setdiff(epi_vars, c(test_accumulators, non_expanded_states))
                 block_size <- length(epi_vars)
               }
             if (!requireNamespace("latticeExtra")) {
               stop("can't add block lines: please install latticeExtra package")
             }
               block_col <- rep(block_col, length.out=length(block_size))
               for (i in seq_along(block_size)) {
                 ## offset by 0.5 so we are illustrating state values
                 block_pos <- seq(nrow(M)+0.5,0,by=-block_size[i])
                 dd <- list(col=block_col[i],pos=block_pos)
                 p <- (p
                   + latticeExtra::layer(lattice::panel.abline(h=pos,col=col), data=dd)
                   + latticeExtra::layer(lattice::panel.abline(v=pos,col=col), data=dd)
                 )
             } ## loop over block_size
         } ## add_blocks
    } else if (method=="igraph") {
       if (!requireNamespace("igraph")) {
           stop("igraph not available")
       } else {
           g <- igraph::graph_from_adjacency_matrix(M)
           plot(g, layout=igraph::layout_as_tree, ...)
       }
    } else if (method=="diagram") {
        xpos <- c(S=1,E=2,Ia=3,Ip=3,Im=4,Is=4,H=5,ICUs=5,ICUd=5,H2=6,D=7,R=7,X=7)
        ypos <- c(S=1,E=1,Ia=1,Ip=2,Im=2,Is=3,H=3,ICUs=4,ICUd=5,H2=4,D=5,R=1,X=4)
        pos <- cbind(xpos,ypos)/8
        M3 <- M[names(xpos),names(xpos)] ## reorder ... does that matter?
        if (do_symbols) {
          M3 <- adjust_symbols(M3)
        } else {
            M3[M3!="0"] <- ""  ## blank out all labels
        }
        diagram::plotmat(t(M3),pos=pos,name=colnames(M3),box.size=box.size, add=FALSE, ...)
    }
    return(p)
}

##' display flow chart on current graphics device
##' @param params parameters (flowchart is determined by non-zero flows)
##' @param testify add testing flows?
##' @param ageify add age structure?
##' @param method visualization method
##' @param ... parameters to pass to \code{\link{show_ratemat}}
##' @examples
##' vis_model()
##' vis_model(testify=TRUE)
##' vis_model(method="diagram")
##' @export
vis_model <- function(params=read_params("PHAC_testify.csv"), testify=FALSE,
                      ageify=FALSE, method=c("Matrix","diagram","igraph"), ...) {
    method <- match.arg(method)
    ## FIXME: accept method= argument, make const_width = (method=="igraph") ?
    state <- make_state(N=1e6, E0=1, params=params)
    state[] <- 1  ## all population states occupied
    M <- make_ratemat(state,params,do_ICU=TRUE, symbols=(method=="diagram"))
    if (testify) {
        M <- testify(M, params)
    }
    if (ageify) {
        M <- ageify(M)
    }
    p <- show_ratemat(M, method=method, ...)
    return(p)
}

adjust_symbols <- function(M) {
    ## use [] throughout to avoid losing dimnames ...
    ## subscripts: _ + letter at end of line
    M[] <- gsub("_([a-z])$","[\\1]",M)
    ## subscripts (letter + number at end of word or line)
    M[] <- gsub("([a-z])([0-9])(\\W|$)","\\1[\\2]\\3", M)
    ## suppress 'nonhosp_mort'
    M[] <- gsub("\\(1 +- +nonhosp_mort\\) +\\*?","",M)
    M[] <- gsub("nonhosp_mort +\\*?","",M)
    ## M -> X == M -> H
    M[grep("^M",M)] <- M["Is","H"]
    M[grep("beta_vec",M)] <- "sum(beta[j]*I[j])"
    return(M)
}

make_flowchart <- vis_model ## back-compatibility


## identify locations within matrix
## ##' @param value return character (TRUE) or numeric (FALSE) position?
pfun <- function(from, to, mat, value=FALSE, recycle=FALSE) {
    ## <start> + label + (_ or <end>)
    from_pos <- grep(sprintf("^%s(_|$)",from), rownames(mat), value=value)
    to_pos <- grep(sprintf("^%s(_|$)",to), colnames(mat), value=value)
    nf <- length(from_pos)
    nt <- length(to_pos)
    if (recycle && (nt==1 || nf==1)) {
        from_pos <- rep(from_pos,length.out=max(nt,nf))
        to_pos <- rep(to_pos,length.out=max(nt,nf))
    }
    if (! (length(to_pos) == length(from_pos) &&
           length(to_pos)>0 && length(from_pos)>0)) {  ## must be positive
        stop(sprintf("to_pos, from_pos don't match: from_pos=%s, to_pos=%s",
                     paste(colnames(mat)[from_pos],collapse=", "),
                     paste(rownames(mat)[to_pos],collapse=", ")
                     ))
    }

    return(cbind(from_pos, to_pos))
}

## exclude states by regex
exclude_states <- function(nm,exclude_states) {
    x_regex <- sprintf("^(%s)_?",paste(exclude_states,collapse="|"))
    xx <- grep(x_regex,nm,invert=TRUE,value=TRUE)
    return(xx)
}

## wrapper for update() to work on lists (using purrr::update_list)
update.list <- function(list, ...) purrr::update_list(list, ...)

## inspired by purrr, infix pkgs
`%||%` <- function (a, b) {
    if (is.null(a)) b else a
}

## logical vector indicating whether first letters of strings are capitalized
first_letter_cap <- function(x) {
  f <- substr(x,1,1)
  return(toupper(f) == f)
}

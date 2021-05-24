## mapping
## (1) dependence of elements of the rate matrix on parameters
##   (so we know which elements of the rate matrix to update when
##    a parameter changes)
## (2) how to construct elements of the rate matrix from
##    parameters (list of indices of multiplied parameters,
##    binary/boolean that specifies whether the parameter's complement
##   (1-x) is taken before multiplying
afun_ind <- function(from, to, val) {
  pos <- pfun(from, to, M)
  v <- substitute(val)
  pars <- vapply(all.vars(v),
                 FUN=match, table=params,
                 FUN.VALUE=numeric(1))
  ## we just need to identify complementarity
  ## do the Bad Thing for now (i.e., deparse + string operations
  ## (alternative is recursive formula processing, ugh.)
  s <- strsplit(deparse(v), "\\*")[[1]]
  complement <- grepl("1 +-",s)
  ratemat_vals <<- append(ratemat_vals, list(pos, pars, complement))
  for (p in pars) {
    ## not quite sure how to make this work ...
    ## ratemat_deps[[p]] <<- append(ratemate_deps[[p]], pos)
    ## ratemat_deps[[p]] <- append(ratemate_deps[[p]], list(pos))
    ## error in append_ratemat(...): object 'ratemate_deps' not found
    ## some weird scoping thing going on here with the `[[<-` operator?
    ratemate_deps <- ratemat_deps
    ratemat_deps[[p]] <- c(ratemate_deps[[p]], list(pos))
  }
  assign("ratemat_deps", ratemat_deps, parent.frame())
  ## add these to a global list
  return(NULL)
}

params <- c("x","y","z")
M <- matrix(0,5,5,dimnames=list(LETTERS[1:5], LETTERS[1:5]))
pfun <- McMasterPandemic:::pfun

## initialize data structures for dependence
ratemat_deps <- vector("list", length(params))
ratemat_vals <- list()

afun_ind("A", "B", x*y*(1-z))
afun_ind("B", "C", z*(1-y))

## construct dependence matrix for graphs?


## FIXME: think about how expansion and accumulation will work together accross
## sub-models (e.g. ageify, testify, vaccify)

## FIXME: document these for real!
##' global variables for testify expansion
#' @export
non_expanded_states <- c("D","X")

##' @rdname non_expanded_states
##' @export
test_extensions <- c("u","p","n","t")

##' @rdname non_expanded_states
##' @export
test_accumulators <- c("N","P")

##' @rdname non_expanded_states
##' @export
asymp_cat <- c("S","E","Ia","Ip","R")

##' @rdname non_expanded_states
severe_cat <- c("Is","H","H2","ICUs","ICUd")

##' @rdname non_expanded_states
## these are 'asymptomatic' (= pre- or asymptomatic)
##' @export
cryptic_cat <- c("Ia","Ip")

#' ## FIXME: document these for real!
#' ##' global variables for vaccify expansion
#' #' @export
#' non_expanded_states <- c()

##' @rdname non_expanded_states
##' @export
vacc_extensions <- c("unvacc","1dose","2dose")

##' @rdname non_expanded_states
##' @export
vacc_accumulators <- c("1dose","2dose")

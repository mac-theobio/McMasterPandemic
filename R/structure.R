## FIXME: think about how expansion and accumulation will work together accross
## sub-models (e.g. ageify, testify, vaccify)

## FIXME: document these for real!
##' global variables for testify expansion
##' @family classic_macpan
##' @export
non_expanded_states <- c("D", "X", "V")

##' @rdname non_expanded_states
##' @family classic_macpan
##' @export
test_extensions <- c("u", "p", "n", "t")

##' @rdname non_expanded_states
##' @family classic_macpan
##' @export
test_accumulators <- c("N", "P")

##' @rdname non_expanded_states
##' @family classic_macpan
##' @export
asymp_cat <- c("S", "E", "Ia", "Ip", "R")

##' @rdname non_expanded_states
severe_cat <- c("Is", "H", "H2", "ICUs", "ICUd")

##' @rdname non_expanded_states
## these are 'asymptomatic' (= pre- or asymptomatic)
##' @family classic_macpan
##' @export
cryptic_cat <- c("Ia", "Ip")

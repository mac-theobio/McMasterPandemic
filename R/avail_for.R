#' What Variables are Available?
#'
#' Find out what variables are available for particular functions that
#' add variables to a model, based on existing variables.
#'
#' \describe{
#'   \item{\code{avail_for_rate}}{\code{\link{add_rate}},\code{\link{rep_rate}},\code{\link{vec_rate}}}
#'   \item{\code{avail_for_sum}}{\code{\link{add_state_param_sum}}}
#'   \item{\code{avail_for_factr}}{\code{\link{add_factr}},\code{\link{vec_factr}}}
#'   \item{\code{avail_for_expr}}{\code{\link{add_sim_report_expr}}}
#'   \item{\code{avail_for_lag}}{\code{\link{add_lag_diff}}}
#'   \item{\code{avail_for_conv}}{\code{\link{add_conv}}}
#' }
#'
#' @param model a \code{\link{flexmodel}} object
#'
#' @return character vector of available variable names
#'
#' @rdname avail_for
#' @export
avail_for_rate = function(model) {
  c(names(model$state),
    names(model$params),
    names(model$sums),
    names(model$factrs))
}

#' @rdname avail_for
#' @export
avail_for_sum = function(model) {
  c(names(model$state),
    names(model$params))
}

#' @rdname avail_for
#' @export
avail_for_factr = function(model) {
  c(names(model$state),
    names(model$params),
    names(model$sums),
    names(model$factrs))
}

#' @rdname avail_for
#' @export
avail_for_pow = function(model) {
  c(names(model$state),
    names(model$params),
    names(model$sums),
    names(model$factrs),
    names(model$pows))
}

#' @rdname avail_for
#' @export
avail_for_expr = function(model) {
  initial_sim_report_names(model)
}

#' @rdname avail_for
#' @export
avail_for_conv = function(model) {
  intermediate_sim_report_names(model)
}

#' @rdname avail_for
#' @export
avail_for_lag = function(model) {
  intermediate_sim_report_names(model)
}

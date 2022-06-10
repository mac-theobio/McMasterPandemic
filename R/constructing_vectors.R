##' Layered Zero State Vector
##'
##' Initialize a state vector with all zeros that is based on
##' a fully-crossed set of states in several sub-models.
##'
##' @param ... character vectors with sub-model state names
##'
##' @export
layered_zero_state = function(...) {
  state_nms = (list(...)
   %>% lapply(as.character)
   %>% Reduce(f = expand_names)
  )
  const_named_vector(state_nms, 0)
}

##' Constant Named Vector
##'
##' @param nms names of the output vector
##' @param cnst single numeric value to be used in every element of the
##' output vector
##'
##' @export
const_named_vector = function(nms, cnst) {
  stopifnot(length(cnst) == 1L)
  stopifnot(is.character(nms))
  setNames(rep(cnst[[1]], length(nms)), nms)
}

#' Merge One Vector into Another by Name
#'
#' If an item in \code{u} has the same name as an item
#' in \code{v} then replace the value in \code{v} with that
#' in \code{u}, otherwise create a new element in \code{v}.
#'
#' @param v named vector or list
#' @param u named vector of list
#' @export
merge_named_vectors = function(v, u) {
  if (is.null(names(v)) | is.null(names(u)) ) {
    stop("v and u must be named vectors")
  }
  for (nm in names(u)) {
    v[nm] = u[nm]
  }
  return(v)
}

##' Merge Named Vector Attribute
##'
##' @param v vector with an attribute, \code{a}
##' @param u vector with attribute \code{a} that will be merged into
##' \code{a} in \code{v} using \code{\link{merge_named_vectors}}
##' @param a name of the atribute to merge
##'
##' @export
merge_named_vec_attr = function(v, u, a) {
  attr_v = attributes(v)
  attr_u = attributes(u)
  attr_v[[a]] = merge_named_vectors(attr_v[[a]], attr_u[[a]])
  attributes(v) = attr_v
  return(v)
}

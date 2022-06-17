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

##' Expand Strain Names
##'
##' This is similar to \code{\link{expand_names}} but for multi-strain models.
##' The expansion also allows setting constraints on the numbers of particular
##' \code{base_states} in the final constrained set of states.
##'
##' @param n_strains number of strains
##' @param base_states character vector giving the states of the base model
##' @param infected_states character vector giving the infected states
##' @param strain_name_prefix prefix for names of strains (they are numbered)
##' @return \code{expand_strain_frame} returns a data frame with one column
##' for each strain and one row for each expanded state.
##' \code{expand_strain_names} returns a vector with the full state names for
##' each state in the expanded model.
expand_strain_frame = function(
    base_states = c("S", "I", "R"),
    constrained_states = c(),
    constraint_counts = 0:1,
    n_strains = 2,
    strain_name_prefix = "strain"
  ) {

  # this implementation creates the full product model and then
  # substracts states that are not in the full model.
  # this implementation could be too slow when the full and
  # constrained models are of very different sizes.


  stopifnot(length(strain_name_prefix) == 1L)

  prod_states = (base_states
    %>% list
    %>% rep(n_strains)
    %>% setNames(strain_name_prefix %_% seq_len(n_strains))
    %>% c(list(stringsAsFactors = FALSE))
    %>% do.call(what = expand.grid)
  )
  no_constraints_mat = matrix(
    as.matrix(prod_states) %in% constrained_states,
    nrow = nrow(prod_states)
  )
  no_constraints = (no_constraints_mat
    %>% rowSums
    %in% constraint_counts
  )
  prod_states[no_constraints, ]
}

##' @rdname expand_strain_frame
##' @inheritDotParams expand_strain_frame
##' @param sep separator between strains
expand_strain_names = function(..., sep = "") {
  unname(apply(expand_strain_frame(...), 1, paste0, collapse = sep))
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

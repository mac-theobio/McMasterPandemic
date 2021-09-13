# Spec Version Management
#
# Functions for managing different versions of specifications
# for flexible modelling.
#
# Specification Document: https://canmod.net/misc/flex_specs


#' Spec Version
#'
#' Return the specification version being assumed by the package.
#' Equivalent to \code{parse_version(getOption('MP_flex_spec_version'))}.
#'
#' Specification Document: https://canmod.net/misc/flex_specs.
#'
#' \code{spec_ver_eq}, \code{spec_ver_gt}, \code{spec_ver_lt}, and
#' \code{spec_ver_btwn} return logical values indicating whether or not
#' the \code{spec_version} is equal to, greater than, less than or between
#' a user-specified version
#'
#' @param version string with the version of the spec
#' (e.g. \code{"0.0.5"})
#' @return \code{svlist} object with the specification version.
#' @importFrom parse_version semver
#' @export
spec_version = function() {
  # https://canmod.net/misc/flex_specs
  parse_version(getOption('MP_flex_spec_version'))
}

#' @rdname spec_version
#' @export
spec_ver_eq = function(version) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  parse_version(current_version) == parse_version(version)
}

#' @rdname spec_version
#' @export
spec_ver_gt = function(version) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  parse_version(current_version) > parse_version(version)
}

#' @rdname spec_version
#' @export
spec_ver_lt = function(version) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  parse_version(current_version) < parse_version(version)
}

#' @rdname spec_version
#' @export
spec_ver_btwn = function(version_left, version_right) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  spec_ver_gt(version_left) &
    spec_ver_lt(version_right)
}

#' Spec and Feature Checks
#'
#' \code{spec_check} will throw an error if a feature is being used that is
#' not supported yet. \code{feature_check} will throw an error if a feature
#' is not being used, even though it is required by the spec version being used.
#'
#' @param introduced_version string with the version of the spec
#' (e.g. \code{"0.0.5"})
#' @param feature free-form text describing the feature
#' @return No return value. Called to get specific error messages when required.
#' @export
spec_check = function(introduced_version, feature) {

  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")

  # TODO: need to handle possessive case?
  verb = ifelse(grepl("s$", feature, perl = TRUE), " were ", " was ")

  if(parse_version(current_version) < parse_version(introduced_version)) {
    stop(
      "\n\n",
      feature, verb, "not introduced until specification version ",
      introduced_version, ".\n",
      "The specification version currently being used is ",
      getOption('MP_flex_spec_version'), '\n',
      "See ", getOption('MP_flex_spec_doc_site'),
      " for more information on specification versions.",
      call. = FALSE)
  }
}

#' @rdname spec_check
#' @export
feature_check = function(introduced_version, feature) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  if(parse_version(current_version) >= parse_version(introduced_version)) {
    stop(
      "\n\n", feature,
      " must be used for all specification versions greater than or equal to ",
      introduced_version, ".\n",
      "The specification version currently being used is ",
      getOption('MP_flex_spec_version'), '\n',
      "See ", getOption('MP_flex_spec_doc_site'),
      " for more information on specification versions.",
      call. = FALSE)
  }
}

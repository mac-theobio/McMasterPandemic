#' String Concatenation Operators
#'
#' Paste with an underscore separator, tilde, or blank space,
#' except for length-zero character strings in the first vector
#'
#' @param x character vector
#' @param y character vector
#'
#' @rdname string_concatenation_operators
#' @export
`%_%` = function(x, y) {
  x = ifelse(nchar(trimws(as.character(x))) == 0L, '', paste(x, '_', sep = ""))
  paste(x, y, sep = "")
}

#' @rdname string_concatenation_operators
#' @export
`%~%` = function(x, y) {
  paste(x, y, sep = " ~ ")
}

#' @rdname string_concatenation_operators
#' @export
`%+%` = function(x, y) paste(x, y, sep = "")


#' Get Substrings by Indices and Separators
#'
#' For example \code{index_sep('a_b_c', 2, '_')} equals \code{'b'}.
#' For example \code{index_sep('a_b_c', c(1, 3), '_')} equals \code{'a_c'}.
#' For example \code{index_sep('a_b_c', -2, '_')} equals \code{'a_c'}.
#' For example \code{index_sep('a_b_c', 4, '_')} equals \code{''}.
#' For example \code{index_sep('a', 1, '_')} equals \code{'a'}.
#' For example \code{index_sep('a', 2, '_')} equals \code{''}.
#' For example \code{index_sep(c('a_b', 'c'), 2, '_')} equals \code{c('b', '')}.
#' For example \code{index_sep('a_b_c', c(3, 1), '_')} equals \code{'c_a'}.
#'
#' @param x character vector
#' @param i integer vector without sign mixing
#' @param sep length-one character vector
#' @param complement if \code{TRUE} the indices in \code{i} that are not matched are returned
#'
#' @export
index_sep = function(x, i, sep = "_") {
  complement = FALSE
  if (any(i < 0L)) {
    if (!all(i < 0L)) stop("cannot mix positive and negative indices")
    complement = TRUE
    i = -1 * i
  }
  stopifnot(length(sep) == 1L)
  if (complement) {
    n_separated_items = nchar(x) - nchar(gsub(sep, '', x)) + 1
    if (length(n_separated_items) > 1L) {
      stop('cannot use complement method with multiple inputs')
    }
    i = setdiff(seq_len(n_separated_items), i)
  }
  (x
   %>% as.character
   %>% strsplit(sep)
   %>% lapply(function(x) {
     ifelse(
       length(x) == 0L,
       '',
       paste0(omit_empty(x[i]), collapse = sep)
     )
   })
   %>% unlist
   %>% unname
  )
}

names_or_values = function(x) {
  if (!is.character(x)) {
    x = names(x)
  }
  x
}

##' Regular Expression Conveniences
##'
##' Specifications of \code{\link{flexmodel}} objects often
##' involves string manipulations with regular expressions,
##' and these functions provide some convenience in this area.
##'
##' @param x character vector to convert into a length-1
##' character vector representing a regular expression
##' pattern
##' @param exact should the regular expression look for
##' exact matches?
##' @param negate should the regular expression look for
##' stings that do not match?
##' @return a regular expression pattern to be used in functions
##' like \code{\link{grep}} and \code{\link{grepl}}, but also in
##' \code{\link{add_state_param_sum}}
##'
##' @rdname macpan_regex_helpers
##' @export
alt_group = function(x, exact = FALSE, negate = FALSE) {
  x = "(" %+% paste0(x, collapse = "|") %+% ")"
  if (negate) {
    x = "(?!(?:" %+% x %+% ")$).*"
  }
  if (exact) {
    x = "^" %+% x %+% "$"
  }
  x
}

##' @rdname macpan_regex_helpers
##' @export
wrap_exact = function(x) {
  "^" %+% x %+% "$"
}

##' @rdname macpan_regex_helpers
##' @export
all_in = names_or_values

##' @rdname macpan_regex_helpers
##' @export
any_var = function(x) {
  (x
   %>% names_or_values
   %>% alt_group(exact = TRUE)
  )
}

##' @rdname macpan_regex_helpers
##' @export
all_except = function(x) {
  (x
   %>% names_or_values
   %>% alt_group(exact = TRUE, negate = TRUE)
  )
}

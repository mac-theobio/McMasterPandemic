#' Class to Represent Scalar, Vector, or Matrix Structure
#'
#' The simplest way to create such objects is with the \code{\link{struc}}
#' function.
#'
#' @slot v character vector giving the expressions for each element
#' of the matrix structure object
#' @slot dims numeric vector giving the dimensions of the matrix
#' structure object
#' @export
setClass('struc', representation(v = "character", dims = "numeric"))

#' Class to Represent an Expanded Matrix Structure Object
#'
#' @slot l list of struc objects
#' @slot dims numeric vector giving the dimensions of the matrix
#' structure object

setClass('struc_expanded', representation(l = "list", dims = "numeric"), validity = function(object) {
  errors = character()
  if(!all(unlist(lapply(object@l, inherits, 'struc')))) errors = append(errors, "not all elements of l are struc objects")
  if(prod(object@dims) != length(object@l)) errors = append(errors, "dimensions are not consistent with the length of l")
  if(length(errors) == 0L) return(TRUE)
  return(errors)
})

#' Construct a struc Object
#'
#' Create model structure objects (see \code{\link{struc-class}}), which are
#' symbolic scalars, vectors, and matrices with elements that are expressions
#' involving parameters, state variables.
#'
#' \code{scal}, \code{vec}, and \code{mat} are similar to \code{struc}, but
#' ensure that the resulting \code{struc} objects are scalars, vectors, or
#' matrices respectfully.
#'
#' \code{struc} objects can be created in other ways and symbolically
#' manipulated using the functions listed below in the 'see also' section.
#'
#' @param ... character vectors or matrices
#' @return \code{struc} object -- see \code{\link{struc-class}}
#' @family struc_functions
#' @export
struc = function(...) {
  l = list(...)
  if(length(l) == 1L) {
    v = l[[1]]
  } else {
    v = unlist(l)
  }
  stopifnot(is.character(v))  # S4 takes care of this check, right?
  stopifnot(is.matrix(v) | is.atomic(v))
  d = get_dim(v)
  dim(v) = NULL
  new('struc', v = wrap_paren(v), dims = d)
}

#' @rdname struc
#' @export
scal = function(...) {
  x = struc(...)
  stopifnot(is_1by1(x))
  x
}

#' @rdname struc
#' @export
vec = struc

#' @rdname struc
#' @export
mat = function(...) {
  x = list(...)
  stopifnot(length(unique(sapply(x, length))) == 1L)
  nrow = length(x)
  ncol = length(x[[1]])
  x = unlist(lapply(x, as.character))
  struc(matrix(x, nrow, ncol, byrow = TRUE))
}

#' Cross Two Sets of Names in Matrix Form
#'
#' @param x character vector
#' @param y character vector
#'
#' @family struc_functions
#' @return
#' @export
cross_mat = function(x, y, sep = "_") {
  struc(matrix(
    expand_names(x, y, sep),
    nrow = length(x),
    ncol = length(y)
  ))
}

#' @rdname struc
#' @param x object to convert to a struc object
#' @export
as.struc = function(x) {
  if(inherits(x, 'struc')) return(x)
  return(struc(x)) # still might fail
}

#' Construct a struc_expanded Object
#'
#' @param l list of struc objects
#' @param d dimensions of the resulting object
#' @return struc_expanded object

struc_expanded = function(l, d) {
  stopifnot(is.recursive(l))
  l = lapply(l, as.struc)
  new('struc_expanded', l = l, dims = d)
}

#' Expand and Contract struc and struc_expanded Objects
#'
#' @param x object to contract or expand
#' @rdname expand_contract_struc

expand_struc = function(x) {
  (x
   %>% as.matrix
   %>% strsplit("\\+")
   %>% lapply(trimws)
   %>% struc_expanded(d = dim(x))
  )
}

contract_struc = function(x) {
  (x
   %>% slot('l')
   %>% lapply(sum)
   %>% lapply(as.character)
   %>% unlist
   %>% matrix(nrow = nrow(x), ncol = ncol(x))
   %>% struc
  )
}

# get dimensions of a matrix or vector
get_dim = function(x) {
  if(is.null(dim(x))) {
    return(c(length(x), 1))
  } else {
    stopifnot(length(dim(x)) == 2L)
  }
  dim(x)
}

# Not used anywhere?
num_prod = function(x) {
  sapply(expand_struc(x)@l, nrow)
}

#' Test if 1-by-1
#'
#' Test if a \code{\link{struc-class}} object is a scalar
#'
#' @param x struc object
#' @family struc_info_functions
#' @return TRUE or FALSE
#' @export
is_1by1 = function(x) {
            all(dim(x) == 1L)
}

#' Test if struc Objects have the Same Dimensions
#'
#' @param x,y \code{\link{struc}} objects
#' @return TRUE or FALSE
#' @family struc_info_functions
#' @export
same_dims = function(x, y) {
            all(dim(x) == dim(y))
          }


setGeneric('resolve',
           function(x) {
             standardGeneric('resolve')
           })


setMethod('resolve', c(x = 'struc'),
          function(x) {
            contract_struc(resolve(expand_struc(x)))
          })


setMethod('resolve', c(x = 'struc_expanded'),
          function(x) {
            (x
             %>% slot('l')
             %>% lapply(gsub, pattern = '(\\(1\\)\\s*\\*{1}|\\*{1}\\s*\\(1\\))', replacement = '')
             %>% lapply(trimws)
             %>% lapply(function(y) ifelse(y == '', '(1)', y))
             %>% lapply(sub, pattern = ".*\\(0\\).*", replacement = "(0)")
             %>% lapply(function(y) {
               y = y[y != "(0)"]
               if(length(y) == 0L) y = "(0)"
               y
             })
             %>% struc_expanded(d = dim(x))
            )
          })

#' Convert struc Object to a matrix
#'
#' @param x \code{\link{struc-class}} object
#' @param ... for S3 method consistency
#' @return matrix
#' @export
as.matrix.struc = function(x, ...) {
  matrix(x@v, x@dims[1], x@dims[2])
}

#' Convert struc Object to a character vector
#'
#' @param x \code{\link{struc-class}} object
#' @param ... for S3 method consistency
#' @return character vector
#' @export
as.character.struc = function(x, ...) {
  as.character(as.matrix(x))
}

pluralizer = function(number, noun) {
  number = as.integer(round(number))
  paste(as.character(number), ' ', noun,
        ifelse(number == 1L, '', 's'), sep = '')
}

wrap_paren = function(x) {
  # not perfect but provides some convenience
  i = which_unwraped(x)
  x[i] = paste0('(', x[i], ')')
  x
}

which_unwraped = function(x) {
  # test for no parentheses, + or * signs
  #  -- this is a weak test, but ok for now
  !grepl('(\\(|\\)|\\+|\\*)', x)
}

#' Symbolic Complement and Inverse
#'
#' @param x character vector
#' @export
#' @examples
#' complement('p')
#' inverse('N')
complement = function(x) {
  i = which_unwraped(x)
  x[i] = "1 - " %+% x[i]
  x
}

#' @rdname complement
#' @export
inverse = function(x) {
  i = which_unwraped(x)
  x[i] = "1 / " %+% x[i]
  x
}


setMethod(f = "show",
          signature = "struc",
          definition = function(object){
            is_one_row = nrow(object) == 1L
            is_one_col = ncol(object) == 1L
            cat(
              "struc object with ",
              pluralizer(nrow(object), "row"),
              " and ",
              pluralizer(ncol(object), "column"),
              ":\n\n", sep = '')
            for(i in 1:nrow(object)) {
              for(j in 1:ncol(object)) {
                if(!is_one_row) cat('Row', i, '\n')
                if(!is_one_col) cat('Column', j, '\n')
                cat(trimws(gsub('\\+', '\\+\n', as.matrix(object)[i, j])),
                    '\n', sep = '')
              }
            }
          })

#' @export
setMethod("*", c(e1 = 'struc_expanded', e2 = 'struc'),
          function(e1, e2) {
            l = lapply(e1@l, simple_mult, e2)
            struc_expanded(lapply(l, struc), d = dim(e1))
          })

#' @export
setMethod('*', c(e1 = 'struc_expanded', e2 = 'struc_expanded'),
          function(e1, e2) {
            # TODO: check that dims match
            l = mapply(simple_mult, e1@l, e2@l, SIMPLIFY = FALSE)
            struc_expanded(lapply(l, struc), d = dim(e1))
          })

#' @describeIn struc Elementwise or scalar multiplication
#' @export
setMethod("*", c(e1 = 'struc', e2 = 'struc'),
          function(e1, e2) {
            if(!same_dims(e1, e2)) {
              if(is_1by1(e1)) {
                big = expand_struc(e2)
                # small always gets expanded to fix bug below
                small = expand_struc(e1)
              } else if(is_1by1(e2)) {
                big = expand_struc(e1)
                # small always gets expanded to fix bug below
                small = expand_struc(e2)
              } else {
                stop('if dims are different, one operand needs to be 1-by-1')
              }
            } else {
              big = expand_struc(e1)
              small = expand_struc(e2)
            }
            # FIXED: this used to fail silently if small is 1-by-1 with more than
            #        one product.
            #        work around was to use kronecker instead of * ...
            #        not sure why this worked
            contract_struc(resolve(big * small))
          }
)

simple_mult = function(x, y) {
  c(outer(x@v, y@v, paste, sep = ' * '))
}

#' @describeIn struc Elementwise or scalar addition
#' @export
setMethod("+", c(e1 = 'struc', e2 = 'struc'),
          function(e1, e2) {
            if(!same_dims(e1, e2)) {
              stopifnot(is_1by1(e1) | is_1by1(e2))
            }
            struc(paste(e1@v, e2@v, sep = ' + '))
          }
)

#' @describeIn struc Matrix multiplication
#' @export
setMethod("%*%", c(x = 'struc', y = 'struc'),
          function(x, y) {
            stopifnot(x@dims[2] == y@dims[1])
            r = x@dims[1]
            c = y@dims[2]
            x = as.matrix(x)
            y = as.matrix(y)
            result = matrix('', r, c)
            for(i in 1:r) {
              for(j in 1:c) {
                xi = struc(x[i, ])
                yj = struc(y[, j])
                result[i, j] = sum(xi * yj)@v
              }
            }
            struc(result)
          }
)

#' @describeIn struc Kronecker product
#' @export
setMethod("kronecker", c(X = 'struc', Y = 'struc'),
          function(X, Y) {
            struc_stretch(X, nrow(Y), ncol(Y)) * struc_block(Y, nrow(X), ncol(X))
          })

#' @describeIn struc Matrix or vector transpose
#' @export
setMethod("t", c(x = 'struc'),
          function(x) {
            x@v = c(t(as.matrix(x)))
            x@dims = rev(x@dims)
            x
          })

setGeneric('inner',
           function(x, y, ...) {
             standardGeneric('inner')
           })

#' @describeIn struc Inner product of vectors or matrices
#' @export
setMethod("inner", c(x = 'struc', y = 'struc'),
    function(x, y) {
      stopifnot(same_dims(x, y))
      sum(x * y)
})

#' @describeIn struc Sum of vector or matrix elements
#' @export
setMethod("sum", c(x = 'struc'),
          function(x, ..., na.rm = FALSE) {
            l = unlist(lapply(c(list(x), list(...)), slot, 'v'))
            struc(paste(l, collapse = ' + '))
          })

#' @describeIn struc Product of vector or matrix elements
#' @export
setMethod("prod", c(x = 'struc'),
          # FIXME: does this work right?
          function(x, ..., na.rm = FALSE) {
            l = unlist(lapply(c(list(x), list(...)), slot, 'v'))
            struc(paste(l, collapse = ' * '))
          })

#' @describeIn struc Row sums of matrices
#' @export
setMethod("rowSums", c(x = 'struc'),
          function(x, na.rm = FALSE, dims = 1L) {
            struc(apply(as.matrix(x), 1, function(y) sum(struc(y))@v))
          })

#' @describeIn struc Column sums of matrices
#' @export
setMethod("colSums", c(x = 'struc'),
          function(x, na.rm = FALSE, dims = 1L) {
            t(struc(apply(as.matrix(x), 2, function(y) sum(struc(y))@v)))
          })

#' @describeIn struc Dimensions of a matrix
#' @export
setMethod("dim", c(x = 'struc'),
          function(x) {
            x@dims
          })


# @describeIn struc Number of matrix rows
# @export
#nrow.struc <- function(x) {
#  x@dims[1]
#}


# @describeIn struc Number of matrix columns
# @export
#setMethod("ncol", c(x = 'struc'),
#          function(x) {
#            x@dims[2]
#          })


setMethod("dim", c(x = 'struc_expanded'),
          function(x) {
            x@dims
          })


# setMethod("nrow", c(x = 'struc_expanded'),
#           function(x) {
#             x@dims[1]
#           })


# setMethod("ncol", c(x = 'struc_expanded'),
#           function(x) {
#             x@dims[2]
#           })

#' Repeat a Block
#'
#' @param x \code{\link{struc-class}} object representing the block
#' @param row_times number of times to replicate the block downwards
#' @param col_times number of times to replicate the block to the right
#' @return struc object
#' @family struc_functions
#' @export
struc_block = function(x, row_times, col_times) {
  stopifnot(is(x, 'struc'))
  x = as.matrix(x)
  x =   matrix(rep(c(  x ), times = col_times), nrow(x), col_times * ncol(x))
  x = t(matrix(rep(c(t(x)), times = row_times), ncol(x), row_times * nrow(x)))
  return(struc(x))
}

#' Stretch a Block
#'
#' @param x struc object representing the block
#' @param row_each number of times to replicate each element of the block downwards
#' @param col_each number of times to replicate each element of the block to the right
#' @return struc object
#' @family struc_functions
#' @export
struc_stretch = function(x, row_each, col_each) {
  stopifnot(is(x, 'struc'))
  x = as.matrix(x)
  x =   matrix(rep(c(  x ), each = row_each), row_each * nrow(x), ncol(x))
  x = t(matrix(rep(c(t(x)), each = col_each), col_each * ncol(x), nrow(x)))
  return(struc(x))
}

#' Numerically Evaluate Struc Object
#'
#' @param x \code{\link{struc-class}} object
#' @param values object coercible to a named list
#'
#' @export
struc_eval = function(x, values) {
  x = as.matrix(x)
  y = matrix(nrow = nrow(x), ncol = ncol(x))
  for(i in 1:nrow(x)) {
    for(j in 1:ncol(x)) {
      y[i, j] = with(as.list(values), eval(parse(text = x[i, j])))
    }
  }
  y
}

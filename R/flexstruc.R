#' Class to Represent Matrix Structure
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
#' @export
setClass('struc_expanded', representation(l = "list", dims = "numeric"), validity = function(object) {
  errors = character()
  if(!all(unlist(lapply(object@l, inherits, 'struc')))) errors = append(errors, "not all elements of l are struc objects")
  if(prod(object@dims) != length(object@l)) errors = append(errors, "dimensions are not consistent with the length of l")
  if(length(errors) == 0L) return(TRUE)
  return(errors)
})

#' Construct a struc Object
#'
#' @param ... character vectors or matrices
#' @return struc object
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
#' @export
struc_expanded = function(l, d) {
  stopifnot(is.recursive(l))
  l = lapply(l, as.struc)
  new('struc_expanded', l = l, dims = d)
}

#' Expand and Contract struc and struc_expanded Objects
#'
#' @param x object to contract or expand
#' @rdname expand_contract_struc
#' @export
expand_struc = function(x) {
  (x
   %>% as.matrix
   %>% strsplit("\\+")
   %>% lapply(trimws)
   %>% struc_expanded(d = dim(x))
  )
}

#' @rdname expand_contract_struc
#' @export
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

#' Test if 1-by-1
#'
#' @param x struc object
#' @return TRUE or FALSE
#' @export
is_1by1 = function(x) {
            all(dim(x) == 1L)
          }

#' Test if struc Objects have the Same Dimensions
#'
#' @param x,y struc objects
#' @return TRUE or FALSE
#' @export
same_dims = function(x, y) {
            all(dim(x) == dim(y))
          }


setGeneric('resolve',
           function(x) {
             standardGeneric('resolve')
           })

#' @export
setMethod('resolve', c(x = 'struc'),
          function(x) {
            contract_struc(resolve(expand_struc(x)))
          })

#' @export
setMethod('resolve', c(x = 'struc_expanded'),
          function(x) {
            (x
             %>% slot('l')
             %>% lapply(gsub, pattern = '\\**\\s*\\(1\\)\\s*\\**', replacement = '')
             %>% lapply(trimws)
             %>% lapply(function(y) ifelse(y == '', '(1)', y))
             %>% struc_expanded(d = dim(x))
            )
          })

#' Convert struc Object to a matrix
#'
#' @param x struc object
#' @return matrix
#' @export
as.matrix.struc = function(x) {
  matrix(x@v, x@dims[1], x@dims[2])
}

#' Convert struc Object to a character vector
#'
#' @param x struc object
#' @return character vector
#' @export
as.character.struc = function(x) {
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

setMethod(f = "show",
          signature = "struc",
          definition = function(object){
            cat(
              "struc object with ",
              pluralizer(nrow(object), "row"),
              " and ",
              pluralizer(ncol(object), "column"),
              ":\n\n", sep = '')
            for(i in 1:nrow(object)) {
              for(j in 1:ncol(object)) {
                cat('Row', i, '\n')
                cat('Column', j, '\n')
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

#' @export
setMethod("*", c(e1 = 'struc', e2 = 'struc'),
          function(e1, e2) {
            if(!same_dims(e1, e2)) {
              if(is_1by1(e1)) {
                big = expand_struc(e2)
                small = e1
              } else if(is_1by1(e2)) {
                big = expand_struc(e1)
                small = e2
              } else {
                stop('if dims are different, one operand needs to be 1-by-1')
              }
            } else {
              big = expand_struc(e1)
              small = expand_struc(e2)
            }
            contract_struc(resolve(big * small))
          }
)

simple_mult = function(x, y) {
  c(outer(x@v, y@v, paste, sep = ' * '))
}

#' @export
setMethod("+", c(e1 = 'struc', e2 = 'struc'),
          function(e1, e2) {
            if(!same_dims(e1, e2)) {
              stopifnot(is_1by1(e1) | is_1by1(e2))
            }
            struc(paste(e1@v, e2@v, sep = ' + '))
          }
)

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

#' @export
setMethod("kronecker", c(X = 'struc', Y = 'struc'),
          function(X, Y) {
            struc_stretch(X, nrow(Y), ncol(Y)) * struc_block(Y, nrow(X), ncol(X))
          })

#' @export
setMethod("t", c(x = 'struc'),
          function(x) {
            x@v = c(t(as.matrix(x)))
            x@dims = rev(x@dims)
            x
          })

#' @export
setMethod("sum", c(x = 'struc'),
          function(x, ..., na.rm = FALSE) {
            l = unlist(lapply(c(list(x), list(...)), slot, 'v'))
            struc(paste(l, collapse = ' + '))
          })

#' @export
setMethod("prod", c(x = 'struc'),
          # FIXME: does this work right?
          function(x, ..., na.rm = FALSE) {
            l = unlist(lapply(c(list(x), list(...)), slot, 'v'))
            struc(paste(l, collapse = ' * '))
          })

#' @export
setMethod("rowSums", c(x = 'struc'),
          function(x, na.rm = FALSE, dims = 1L) {
            struc(apply(as.matrix(x), 1, function(y) sum(struc(y))@v))
          })

#' @export
setMethod("colSums", c(x = 'struc'),
          function(x, na.rm = FALSE, dims = 1L) {
            t(struc(apply(as.matrix(x), 2, function(y) sum(struc(y))@v)))
          })

#' @export
setMethod("dim", c(x = 'struc'),
          function(x) {
            x@dims
          })

#' @export
setMethod("nrow", c(x = 'struc'),
          function(x) {
            x@dims[1]
          })

#' @export
setMethod("ncol", c(x = 'struc'),
          function(x) {
            x@dims[2]
          })

#' @export
setMethod("dim", c(x = 'struc_expanded'),
          function(x) {
            x@dims
          })

#' @export
setMethod("nrow", c(x = 'struc_expanded'),
          function(x) {
            x@dims[1]
          })

#' @export
setMethod("ncol", c(x = 'struc_expanded'),
          function(x) {
            x@dims[2]
          })

#' Repeat a Block
#'
#' @param x struc object representing the block
#' @param row_times number of times to replicate the block downwards
#' @param col_times number of times to replicate the block to the right
#' @return struc object
#' @export
struc_block = function(x, row_times, col_times) {
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
#' @export
struc_stretch = function(x, row_each, col_each) {
  x = as.matrix(x)
  x =   matrix(rep(c(  x ), each = row_each), row_each * nrow(x), ncol(x))
  x = t(matrix(rep(c(t(x)), each = col_each), col_each * ncol(x), nrow(x)))
  return(struc(x))
}

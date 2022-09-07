# classes ------------------------------------------------

#' Class to Represent Scalar, Vector, or Matrix Structure
#'
#' The simplest way to create such objects is with the \code{\link{struc}}
#' function.
#'
#' @param e1 first argument of a binary operator
#' @param e2 second argument of a binary operator
#' @param x argument in a method
#' @param y argument in a method
#' @param X argument in a method
#' @param Y argument in a method
#' @param ... arguments to pass on to other methods
#' @param na.rm should missing values be removed?
#' @param dims object dimensions
#'
#' @slot v character vector giving the expressions for each element
#' of the matrix structure object
#' @slot dims numeric vector giving the dimensions of the matrix
#' structure object
#'
#' @note Methods that involve objects of class \code{struc_expanded}
#' in their signature are not intended for users.
#'
#' @export
setClass("struc", representation(v = "character", dims = "numeric"))

setClass(
  "struc_dimnames",
  representation(dimnames = "list"),
  contains = "struc",
  validity = function(object) {
    if (length(object@dimnames[[1]]) != object@dims[1]) {
      return("not the right amount of row names")
    }
    if (length(object@dimnames[[2]]) != object@dims[2]) {
      return("not the right amount of column names")
    }
    if (!all(vapply(object@dimnames, is.character, logical(1L)))) {
      return("dimnames must be character vectors")
    }
    TRUE
  }
)

# Class to Represent an Expanded Matrix Structure Object
#
# @slot l list of struc objects
# @slot dims numeric vector giving the dimensions of the matrix
# structure object

setClass(
  'struc_expanded',
  representation(l = "list", dims = "numeric"),
  validity = function(object) {
  errors = character()
  if (!all(unlist(lapply(object@l, inherits, 'struc')))) errors = append(errors, "not all elements of l are struc objects")
  if (prod(object@dims) != length(object@l)) errors = append(errors, "dimensions are not consistent with the length of l")
  if (length(errors) == 0L) return(TRUE)
  return(errors)
}
)

#' @export
setClass(
  "EpiMatrix",
  representation(
    name = "character",
    default = "matrix",
    elnames = "character"
  ),
  validity = function(object) {
    if (length(object@name) != 1L) {
      return("the name of an epidemiological matrix must be a length-1 character vector")
    }
    if (!good_names(object@name)) {
      return("names of epidemiological matrices must start with a letter and contain only letters, numbers, and underscores")
    }
    if (is.null(dim(object@default))) {
      return("epidemiological matrices must have dimensions")
    }
    if (length(dim(object)) != 2L) {
      return("epidemiological matrices must have two dimensions (i.e. a matrix)")
    }
    if (length(rownames(object)) != nrow(object)) {
      return("number of row names must match the number of rows")
    }
    if (any(is.na(rownames(object)))) {
      return("epidemiological matrices cannot have missing row names")
    }
    if (any(duplicated(rownames(object)))) {
      return("epidemiological matrices cannot have duplicate row names")
    }
    if (length(colnames(object)) != ncol(object)) {
      return("epidemiological matrices number of column names must match the number of columns")
    }
    if (any(is.na(colnames(object)))) {
      return("epidemiological matrices cannot have missing column names")
    }
    if (any(duplicated(colnames(object)))) {
      return("epidemiological matrices cannot have duplicate column names")
    }
    if (length(names(object)) != size(object)) {
      return("number of element names must match the number of elements")
    }
    if (any(is.na(names(object)))) {
      return("epidemiological matrices cannot have missing element names")
    }
    if (any(duplicated(names(object)))) {
      return("epidemiological matrices cannot have duplicate element names")
    }
    if (any(is.na(object@default))) {
      return("default epidemiological matrices cannot have missing values")
    }
  }
)

#' @export
setClass(
  "EpiMatrixInput",
  contains = "EpiMatrix",
  validity = function(object) {
    if (!is.numeric(object@default)) {
      return("input epidemiological matrices must be numeric-valued")
    }
    TRUE
  }
)

#' @export
setClass(
  "EpiMatrixDerived",
  contains = "EpiMatrix",
  validity = function(object) {
    if (!is.character(object@default)) {
      return("derived epidemiological matrices must be character-valued")
    }
    TRUE
  }
)

#' @export
setClass(
  "EpiModelMatrices",
  representation(l = "list"),
  validity = function(object) {
    TRUE
  }
)

# constructors ------------------------------------------------

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
  if (length(l) == 1L) {
    v = l[[1]]
    if (is(v, "EpiMatrixDerived")) {
      v = v@default
    }
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
#' @param sep character for separating \code{x} and \code{y} components
#'
#' @family struc_functions
#' @export
cross_mat = function(x, y, sep = "_") {
  struc(matrix(
    expand_names(x, y, sep),
    nrow = length(x),
    ncol = length(y)
  ))
}

#' @export
epi_mat = function(nm, x) {
  x = as.matrix(x)
  if (length(nm) == 1L)
  if (is.null(rownames(x))) {
    if (nrow(x) != 1L) {
      stop("row names must be supplied for variables with more than one row")
    }
    rownames(x) = nm
  }
  if (is.null(colnames(x))) {
    if (ncol(x) != 1L) {
      stop("column names must be supplied for variables with more than one column")
    }
    colnames(x) = nm
  }

  elnames = make_elnames(x, nm)

  if ((nrow(x) == 1L) & any(rownames(x) != nm)) rownames(x) = nm
  if ((ncol(x) == 1L) & any(colnames(x) != nm)) colnames(x) = nm

  if (is.character(x)) {
    return(new("EpiMatrixDerived", name = nm, default = x, elnames = elnames))
  } else if (is.numeric(x)) {
    mode(x) = "double"
    return(new("EpiMatrixInput", name = nm, default = x, elnames = elnames))
  }
  stop("invalid input")
}

#' @export
epi_mat_list = function(...) {
  mats = unlist(lapply(list(...), list_if_not_list), recursive = FALSE)
  l = mapply(epi_mat, names(mats), mats, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  #nms = vapply(l, slot, "name", character(1L))
  if (!all(vapply(l, is, logical(1L), class2 = "EpiMatrix"))) {
    return("all matrices must be of class EpiMatrix")
  }
  # slot_nms = vapply(l, slot, character(1L), name = "name", USE.NAMES = FALSE)
  # if (!isTRUE(all.equal(names(l), slot_nms))) {
  #   return("the names of list elements must equal the names of the elements themselves")
  # }
  structure(
    setNames(l, names(mats)),
    class = "EpiMatrixList",
    const = list()
  )
}

#' @export
epi_const_list = function(...) {
  l = unlist(lapply(list(...), list_if_not_list), recursive = FALSE)
  if (!all(vapply(l, is_num_or_chr, logical(1L)))) {
    return("all constants must be numeric")
  }
  structure(l, class = "EpiConstList")
}


#' #' @export
#' epi_model_init = structure(
#'   list(
#'     const_list = epi_const_list(),
#'     mat_list = epi_mat_list(),
#'     param_mat_nms = character(0L),
#'     state_mat_nms = character(0L),
#'     rate_list = epi_rate_list()
#'   ),
#'   class = "EpiModelInit"
#' )


# generic definitions ------------------------------------------------

#' Derive
#'
#' Derive symbolic expressions of existing epidemiological variables, and
#' add them to the list
#'
#' @export
derive = function(mat_list, ...) {
  UseMethod("derive")
}

#' @export
const = function(mat_list, ...) {
  UseMethod("const")
}

#' @export
define_state = function(mat_list, ...) {
  UseMethod("define_state")
}

#' Alter
#'
#' Alter existing epidemiological variables, and add them to the list
#'
#' @export
alter = function(mat_list, ...) {
  UseMethod("alter")
}

#' @export
setGeneric("dimnames_from_template", function(x, template) {
  standardGeneric("dimnames_from_template")
})

#' Names to Values
#'
#' Replace the values of an object with its names.
#'
#' This function only makes sense for objects that have the same shape
#' as their names (e.g. vectors, EpiMatrix ... well that might be it actually)
#'
#' @param x object to fill with its names
#'
#' @export
names_to_values = function(x) {
  UseMethod("names_to_values")
}

#' @export
names_from_values = function(x) {
  UseMethod("names_from_values")
}

#' Symbolic Complement and Inverse
#'
#' @param x character vector, \code{\link{struc}} object, or
#' \code{\link{EpiMatrixInput}} object
#' @export
#' @examples
#' complement('p')
#' inverse('N')
setGeneric("complement", function(x) {
  standardGeneric("complement")
})

#' @rdname complement
#' @export
setGeneric("inverse", function(x) {
  standardGeneric("inverse")
})

#' @export
setGeneric("inner", function(x, y, ...) {
  standardGeneric("inner")
})

#' @export
setGeneric("lin_combin", function(w, l) {
  standardGeneric("lin_combin")
})

#' Symbolic Matrix Diagonal
#'
#' Get the diagonal of a \code{\link{struc}} matrix.
#'
#' @param x \code{\link{struc}} object
#'
#' @export
setGeneric("diagonal", function(x) {
  standardGeneric("diagonal")
})

#' Size
#'
#' Generic definition of object size.
#'
#' @param x object
#'
#' @export
setGeneric("size", function(x) {
  standardGeneric("size")
})

#' Flatten Symbolic Matrix
#'
#' @param x \code{\link{struc}} object
#'
#' @export
setGeneric("flatten", function(x) {
 standardGeneric("flatten")
})

#' @export
setGeneric("values", function(x) standardGeneric("values"))

#' Symbolic Row and Column Multiplication
#'
#' Element-wise multiplication of the rows/columns of a matrix by
#' a vector of the same length as each row/column.
#'
#' @param x \code{\link{struc}} object playing the role of a matrix
#' @param y \code{\link{struc}} object playing the role of a vector
#' @name row_col_mult
NULL

#' @rdname row_col_mult
#' @export
setGeneric("row_mult", function(x, y) {
 standardGeneric("row_mult")
})

#' @rdname row_col_mult
#' @export
setGeneric("col_mult", function(x, y) {
 standardGeneric("col_mult")
})

setGeneric('resolve', function(x) {
 standardGeneric('resolve')
})

# exposed utility methods ------------------------------------------------

#' @describeIn struc Extract specific elements of a struc object
#' @export
`[.struc` = function(x, i, j, ..., drop = FALSE) {
  if (length(list(...)) != 0L) {
    stop("struc objects cannot contain multidimensional arrays")
  }
  if (drop) stop("not allowed to drop dimensions in struc subsetting")
  x_mat = as.matrix(x)
  if (is(x, "struc_dimnames")) {
    # this might not be necessary
    dimnames(x_mat) = dimnames(x)
  }
  y_mat = x_mat[i, j, drop = drop]
  y = struc(y_mat)
  if (is(x, "struc_dimnames")) {
    dimnames(y) = dimnames(y_mat)
  }
  y
}

#' @export
`[.EpiMatrix` = function(x, i, j, ..., drop = FALSE) {
  epi_mat(x@name, x@default[i, j, drop = FALSE])
}

epi_mat_replace_utility = function(x, i, j, value) {
  y = x@default
  y[i, j] = c(value)
  epi_mat(x@name, y)
}

#' @export
`[<-.EpiMatrixInput` = function(x, i, j, ..., value) {
  if (!is.numeric(value)) {
    stop("can only replace values of an epidemiological input matrix with numerical values")
  }
  epi_mat_replace_utility(x, i, j, value)
}

#' @export
`[<-.EpiMatrixDerived` = function(x, i, j, ..., value) {
  if (!is.character(value)) {
    stop("can only replace values of an epidemiological input matrix with numerical values")
  }
  epi_mat_replace_utility(x, i, j, value)
}

#' @export
`[.EpiMatrixList` = function(x, i, j, ..., drop = FALSE) {
  cls = class(x)
  y = unclass(x)[i, drop = FALSE]
  class(y) = cls
  y
}

#' @export
setMethod("[[", c(x = "EpiModelMatrices", i = "character"), function(x, i) {
  x@l[[i]]
})

#' @export
as.double.EpiMatrixList = function(x) {
  if (length(select_derived_mats(x)) != 0L) {
    warning("list contained derived matrices, which cannot be converted to numeric")
  }
  x = select_input_mats(x)
  setNames(
    unlist(lapply(x, values), use.names = FALSE, recursive = FALSE),
    unlist(lapply(x, names), use.names = FALSE, recursive = FALSE)
  )
}

#' @describeIn struc Dimensions of a matrix
#' @export
setMethod("dim", c(x = 'struc'), function(x) {
  x@dims
})

#' @export
setMethod("dim", c(x = "EpiMatrix"), function(x) dim(x@default))


#' @export
setMethod("dimnames", c(x = "struc_dimnames"), function(x) {
  x@dimnames
})

#' @export
setMethod("names", c(x = "EpiMatrix"), function(x) x@elnames)

#' @export
setMethod("dimnames<-", c(x = "struc", value = "list"), function(x, value) {
  new("struc_dimnames", v = x@v, dims = x@dims, dimnames = value)
})

#' @export
setMethod("dimnames", c(x = "EpiMatrix"), function(x) dimnames(x@default))

#' @export
setMethod("dimnames_from_template", c(x = "struc", template = "struc_dimnames"), function(x, template) {
  dimnames(x) = dimnames(template)
  x
})

#' @export
names_to_values.numeric = function(x) {
  assert_good_names(x)
  x[] = names(x)
  x
}

#' @export
names_to_values.character = function(x) {
  assert_good_names(x)
  x[] = names(x)
  x
}

#' @export
names_to_values.matrix = function(x) {
  x[] = make_elnames(x, "")
  x
}

#' @export
names_to_values.EpiMatrixDerived = function(x) {
  x[] = names(x)
  x
}

#' @export
names_to_values.EpiMatrixInput = function(x) {
  y = x@default
  y[] = names(x)
  epi_mat(x@name, y)
}

#' @export
names_to_values.array = function(x) stop("this feature is not developed")

#' @export
names_from_values.EpiMatrixDerived = function(x) {
  # could fail, but failure is probably correct as long as the error message
  # is cleaned up. the reason that failure is correct is that not all symbolic
  # matrices have one single name in each cell
  unwrap_paren(as.character(x))
}

#' @export
names_from_values.struc = function(x) {
  # could fail, but that's probably correct as long as the error message
  # is cleaned up
  unwrap_paren(as.character(x))
}

#' @describeIn struc number of elements of \code{struc} object
#' @export
setMethod("size", c(x = "struc"), function(x) {
  length(x@v)
})

#' @export
setMethod("size", c(x = "EpiMatrix"), function(x) prod(dim(x@default)))

#' @export
setMethod("values", c(x = "EpiMatrix"), function(x) unname(flatten(x)))

#' @describeIn struc diagonal of a struc matrix
#' @export
setMethod("diagonal", c(x = "struc"), function(x) {
  struc(diag(as.matrix(x)))
})

#' @describeIn struc Flatten a struc matrix
#' @export
setMethod("flatten", c(x = 'struc'), function(x) {
  struc(as.character(x))
})

#' @export
setMethod("flatten", c(x = "EpiMatrix"), function(x) {
  xx = x@default
  dim(xx) = NULL
  setNames(xx, names(x))
})

#' @rdname complement
#' @export
setMethod("complement", c(x = "character"), function(x) {
  i = which_unwraped(x)
  # TODO: warn if any are wrapped
  x[i] = "1 - " %+% x[i]
  x
})

#' @rdname complement
#' @export
setMethod("complement", c(x = "struc"), function(x) {
  z = complement(unwrap_paren(x@v))
  dim(z) = dim(x)
  struc(z)
})

#' @rdname complement
#' @export
setMethod("complement", c(x = "struc_dimnames"), function(x) {
  y = callNextMethod()
  dimnames(y) = dimnames(x)
  y
})

#' @rdname complement
#' @export
setMethod("complement", c(x = "EpiMatrixInput"), function(x) {
  y = as.matrix(x)
  y[] = complement(names(x))
  epi_mat(x@name %_% "complement", y)
})

#' @rdname complement
#' @export
setMethod("complement", c(x = "numeric"), function(x) {
  if (is.null(names(x))) {
    stop("can only take the symbolic complement of numeric vectors if they are named")
  }
  if (!good_names(names(x))) {
    stop("names are not appropriate for taking symbolic complement")
  }
  complement(names(x))
})

#' @rdname complement
#' @export
setMethod("inverse", c(x = "character"), function(x) {
  i = which_unwraped(x)
  x[i] = "1 / " %+% x[i]
  x
})

#' @rdname complement
#' @export
setMethod("inverse", c(x = "struc"), function(x) {
  z = inverse(unwrap_paren(x@v))
  dim(z) = dim(x)
  struc(z)
})

#' @rdname complement
#' @export
setMethod("inverse", c(x = "struc_dimnames"), function(x) {
  y = callNextMethod()
  dimnames(y) = dimnames(x)
  y
})

#' @rdname complement
#' @export
setMethod("inverse", c(x = "EpiMatrixInput"), function(x) {
  y = as.matrix(x)
  y[] = inverse(names(x))
  epi_mat(x@name %_% "inverse", y)
})

#' @rdname complement
#' @export
setMethod("inverse", c(x = "numeric"), function(x) {
  if (is.null(names(x))) {
    stop("can only take the symbolic inverse of numeric vectors if they are named")
  }
  if (!good_names(names(x))) {
    stop("names are not appropriate for taking symbolic inverse")
  }
  inverse(names(x))
})


setMethod(f = "show", signature = "struc", definition = function(object){
  is_one_row = nrow(object) == 1L
  is_one_col = ncol(object) == 1L
  cat(
    "struc object with ",
    pluralizer(nrow(object), "row"),
    " and ",
    pluralizer(ncol(object), "column"),
    ":\n\n", sep = '')
  for (i in 1:nrow(object)) {
    for (j in 1:ncol(object)) {
      if (!is_one_row) cat('Row', i, '\n')
      if (!is_one_col) cat('Column', j, '\n')
      cat(trimws(gsub('\\+', '\\+\n', as.matrix(object)[i, j])),
          '\n', sep = '')
    }
  }
})

#' @export
setMethod("show", c(object = "EpiMatrix"), function(object) print(object@default))

# exposed utility functions --------------------------

#' @export
select_derived_mats = function(mat_list) {
  mat_list[vapply(mat_list, is, logical(1L), class2 = "EpiMatrixDerived")]
}

#' @export
select_input_mats = function(mat_list) {
  mat_list[vapply(mat_list, is, logical(1L), class2 = "EpiMatrixInput")]
}

#' @export
mats_to_strucs = function(mat_list) {
  struc_list = (mat_list
    %>% lapply(names_to_values)
    %>% lapply(as.struc)
  )
  for (i in seq_along(struc_list)) {
    dimnames(struc_list[[i]]) = dimnames(mat_list[[i]])
  }
  struc_list
}

#' @export
add_derived_matrix = function(model, mat) {
  vec_factr(model, names(mat), values(mat))
}

#' @export
add_derived_matrices = function(model, mat_list) {
  mat_list = select_derived_mats(mat_list)
  if (length(mat_list) == 0L) {
    stop("no derived matrices in mat_list")
  }
  for (mat in mat_list) {
    model = vec_factr(model, names(mat), values(mat))
  }
  model
}

# coercion ------------------------------------------------

#' @rdname struc
#' @param x object to convert to a struc object
#' @export
as.struc = function(x) {
  if (inherits(x, 'struc')) return(x)
  return(struc(x)) # still might fail
}

#' Convert struc Object to a matrix
#'
#' @param x \code{\link{struc-class}} object
#' @param ... for S3 method consistency
#' @return matrix
#' @export
as.matrix.struc = function(x, ...) {
  matrix(x@v, x@dims[1], x@dims[2])
}

#' @export
as.matrix.struc_dimnames = function(x, ...) {
  dn = dimnames(x)
  x = matrix(x@v, x@dims[1], x@dims[2])
  dimnames(x) = dn
  x
}

#' @export
as.matrix.EpiMatrix = function(x) as(x, "matrix")

#' Convert struc Object to a character vector
#'
#' @param x \code{\link{struc-class}} object
#' @param ... for S3 method consistency
#' @return character vector
#' @export
as.character.struc = function(x, ...) {
  as.character(as.matrix(x))
}

#' @export
as.character.EpiMatrix = function(x) as(x, "character")


#' @export
setAs(from = "EpiMatrixDerived", to = "character", function(from) {
  x = as.matrix(from)
  dim(x) = NULL
  x
})

#' @export
setAs(from = "EpiMatrixInput", to = "character", function(from) {
  names(from)
})

#' @export
setAs(from = "EpiMatrix", to = "matrix", function(from) {
  from@default
})

#' @export
setAs(from = "EpiMatrixDerived", to = "struc", function(from) {
  struc(as.character(from))
})

#' @export
setAs(from = "EpiMatrixInput", to = "EpiMatrixDerived", function(from) {
  names_to_values(from)
})

# testing functions ------------------------------------------------

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

#' @rdname is_1by1
#' @export
is_1byn = function(x) {
  (nrow(x) == 1L) & (ncol(x) != 1L)
}

#' @rdname is_1by1
#' @export
is_nby1 = function(x) {
  (nrow(x) != 1L) & (ncol(x) == 1L)
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

#' @rdname same_dims
#' @export
same_dimnames = function(x, y) {
  isTRUE(all.equal(dimnames(x), dimnames(y)))
}

#' @rdname same_dims
#' @export
same_rownames = function(x, y) {
  isTRUE(all.equal(rownames(x), rownames(y)))
}


#' @rdname same_dims
#' @export
same_colnames = function(x, y) {
  isTRUE(all.equal(colnames(x), colnames(y)))
}

# mathematical methods ------------------------------------------------

#' @describeIn struc Elementwise or scalar multiplication
#' @export
setMethod("*", c(e1 = 'struc', e2 = 'struc'), function(e1, e2) {
    if (!same_dims(e1, e2)) {
      if (is_1by1(e1)) {
        big = expand_struc(e2)
        # small always gets expanded to fix bug below
        small = expand_struc(e1)
      } else if (is_1by1(e2)) {
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

#' @export
setMethod("*", c(e1 = "struc_dimnames", e2 = "struc_dimnames"), function(e1, e2) {
  x = callNextMethod()
  dimnames(x) = get_elementwise_dimnames(e1, e2)
  x
})

#' @export
setMethod("+", c(e1 = "struc_dimnames", e2 = "struc_dimnames"), function(e1, e2) {
  x = callNextMethod()
  dimnames(x) = get_elementwise_dimnames(e1, e2)
  x
})

#' @describeIn struc Elementwise or scalar addition
#' @export
setMethod("+", c(e1 = 'struc', e2 = 'struc'),
  function(e1, e2) {
    if (!same_dims(e1, e2)) {
      stopifnot(is_1by1(e1) | is_1by1(e2))
    }
    pp = paste(e1@v, e2@v, sep = ' + ')
    dim(pp) = dim(e1)
    return(struc(pp))
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
            for (i in 1:r) {
              for (j in 1:c) {
                xi = struc(x[i, ])
                yj = struc(y[, j])
                result[i, j] = sum(xi * yj)@v
              }
            }
            struc(result)
          }
)

setMethod("%*%", c(x = "struc_dimnames", y = "struc_dimnames"), function(x, y) {
  z = callNextMethod()
  dimnames(z) = list(rownames(x), colnames(y))
  z
})

#' @describeIn struc Kronecker product
#' @export
setMethod("kronecker", c(X = 'struc', Y = 'struc'),
          function(X, Y) {
            struc_stretch(X, nrow(Y), ncol(Y)) * struc_block(Y, nrow(X), ncol(X))
          })

#' @describeIn struc Matrix or vector transpose
#' @export
setMethod("t", c(x = 'struc'), function(x) {
  x@v = c(t(as.matrix(x)))
  x@dims = rev(x@dims)
  x
})

#' @export
setMethod("t", c(x = "struc_dimnames"), function(x) {
  y = callNextMethod()
  dimnames(y) = rev(dimnames(x))
  y
})

#' @export
setMethod("t", c(x = "EpiMatrix"), function(x) {
  epi_mat(x@name, t(x@default))
})

#' @export
setMethod("inner", c(x = 'struc', y = 'struc'), function(x, y) {
  stopifnot(same_dims(x, y))
  sum(x * y)
})

#' @export
setMethod("lin_combin", c(w = 'struc', l = 'list'), function(w, l) {
  if (!all(vapply(l, is, logical(1L), 'struc'))) {
    stop("linear combinations can only be performed on lists of all struc objects")
  }
  if (length(unique(lapply(l, dim))) != 1L) {
    stop("all elements in the list associated with a linear combination must have the same dimensions")
  }
  if (!all_elementwise_compatible_dimnames(l)) {
    stop("all elements in the list associated with a linear combination must have compatible dimnames")
  }
  if (size(w) != length(l)) {
    stop("the weights, w, must have as many elements as the number of struc objects in the list, l")
  }
  if (length(l) == 1L) return(l[[1L]] * w[1L])
  w = struc_to_list_of_scalars(w)
  wl = mapply(scal_mult, w, l, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  y = wl[[1L]]
  for (i in 2:length(l)) {
    y = y + wl[[i]]
  }
  dimnames(y) = get_elementwise_dimnames(l[[1L]], l[[2L]])
  y
})

#' @describeIn struc Sum of vector or matrix elements
#' @export
setMethod("sum", c(x = 'struc'), function(x, ..., na.rm = FALSE) {
  l = unlist(lapply(c(list(x), list(...)), slot, 'v'))
  struc(paste(l, collapse = ' + '))
})

#' @describeIn struc Product of vector or matrix elements
#' @export
setMethod("prod", c(x = 'struc'), function(x, ..., na.rm = FALSE) {
  # FIXME: does this work right?
  l = unlist(lapply(c(list(x), list(...)), slot, 'v'))
  struc(paste(l, collapse = ' * '))
})

#' @describeIn struc Row sums of matrices
#' @export
setMethod("rowSums", c(x = 'struc'), function(x, na.rm = FALSE, dims = 1L) {
  struc(apply(as.matrix(x), 1, function(y) sum(struc(y))@v))
})

#' @describeIn struc Column sums of matrices
#' @export
setMethod("rowSums", c(x = 'struc_dimnames'), function(x, na.rm = FALSE, dims = 1L) {
  y = callNextMethod()
  dimnames(y) = list(rownames(x), "null")
  y
})

#' @describeIn struc Column sums of matrices
#' @export
setMethod("colSums", c(x = 'struc'), function(x, na.rm = FALSE, dims = 1L) {
  t(struc(apply(as.matrix(x), 2, function(y) sum(struc(y))@v)))
})


#' @describeIn struc Column sums of matrices
#' @export
setMethod("colSums", c(x = 'struc_dimnames'), function(x, na.rm = FALSE, dims = 1L) {
  y = callNextMethod()
  dimnames(y) = list("null", colnames(x))
  y
})


# mathematical functions ------------------------------------------------

#' Repeat a Block
#'
#' @param x \code{\link{struc-class}} object representing the block
#' @param row_times number of times to replicate the block downwards
#' @param col_times number of times to replicate the block to the right
#' @param dimname_template optional object with dimnames to use for the
#' output
#' @return struc object
#' @family struc_functions
#' @export
struc_block = function(x, row_times, col_times, dimname_template = NULL) {
  stopifnot(is(x, 'struc'))
  y = as.matrix(x)
  y =   matrix(rep(c(  y ), times = col_times), nrow(y), col_times * ncol(y))
  y = t(matrix(rep(c(t(y)), times = row_times), ncol(y), row_times * nrow(y)))
  y = struc(y)
  if (!is.null(dimname_template)) {
    y = dimnames_from_template(y, dimname_template)
  }
  return(y)
}

#' Stretch a Block
#'
#' @param x struc object representing the block
#' @param row_each number of times to replicate each element of the block
#' downwards
#' @param col_each number of times to replicate each element of the block
#' to the right
#' @param dimname_template optional object with dimnames to use for the
#' output
#' @return struc object
#' @family struc_functions
#' @export
struc_stretch = function(x, row_each, col_each, dimname_template = NULL) {
  stopifnot(is(x, 'struc'))
  y = as.matrix(x)
  y =   matrix(rep(c(  y ), each = row_each), row_each * nrow(y), ncol(y))
  y = t(matrix(rep(c(t(y)), each = col_each), col_each * ncol(y), nrow(y)))
  y = struc(y)
  if (!is.null(dimname_template)) {
    y = dimnames_from_template(y, dimname_template)
  }
  return(y)
}

#' Block Diagonal Matrix
#'
#' @param x list of \code{struc} objects
#' @return block diagonal struc object
#' @export
struc_bdiag = function(x) {
  y = matrix(
    '(0)',
    sum(unlist(lapply(x, nrow))),
    sum(unlist(lapply(x, ncol)))
  )
  row_pointer = 0L
  col_pointer = 0L
  for (i in seq_along(x)) {
    ii = row_pointer + seq_len(nrow(x[[i]]))
    jj = col_pointer + seq_len(ncol(x[[i]]))
    y[ii, jj] = as.matrix(x[[i]])
    row_pointer = max(ii)
    col_pointer = max(jj)
  }
  return(struc(y))
}

#' @describeIn struc Element-wise multiplication of the rows of a
#' matrix by a vector of the same length as each row.
#' @export
setMethod("row_mult", c(x = "struc", y = "struc"), function(x, y) {
  if (ncol(x) != size(y)) stop("y must be the same size as a row of x")
  z = struc_stretch(flatten(y), nrow(x), 1L) * flatten(x)
  z_values = z@v
  dim(z_values) = dim(x)
  struc(z_values)
})

#' @describeIn struc Element-wise multiplication of the columns of a
#' matrix by a vector of the same length as each column.
#' @export
setMethod("col_mult", c(x = "struc", y = "struc"), function(x, y) {
  if (nrow(x) != size(y)) stop("y must be the same size as a column of x")
  z = struc_block(flatten(y), ncol(x), 1L) * flatten(x)
  z_values = z@v
  dim(z_values) = dim(x)
  struc(z_values)
})



# update epidemiological model -------------------------




# evaluation ----------------------------

#' Numerically Evaluate Struc Object
#'
#' @param x \code{\link{struc-class}} object
#' @param values object coercible to a named list of numeric objects
#'
#' @export
struc_eval = function(x, values) {
  x = as.matrix(x)
  y = matrix(nrow = nrow(x), ncol = ncol(x))
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      y[i, j] = with(as.list(values), eval(parse(text = x[i, j])))
    }
  }
  y
}

#' @export
derive.EpiMatrixList = function(mat_list, ...) {
  valid_funcs = getOption("MP_valid_derive_funcs")
  struc_list = mats_to_strucs(mat_list)
  valid_const = attr(mat_list, "const")

  new_mats = list(...)
  new_mat_names = list()
  for (i in seq_along(new_mats)) {
    mat_i = eval(
      get_rhs(new_mats[[i]]),
      c(struc_list, valid_funcs, valid_const),
      enclos = emptyenv()
    )
    new_mat_names[[i]] = as.character(get_lhs(new_mats[[i]]))
    new_mats[[i]] = epi_mat(
      nm = new_mat_names[[i]],
      x = as.matrix(mat_i)
    )
  }
  names(new_mats) = as.character(new_mat_names)
  structure(
    c(mat_list, new_mats),
    class = "EpiMatrixList",
    const = valid_const
  )
}

#' @export
const.EpiMatrixList = function(mat_list, ...) {
  new_constants = list(...)
  if (!all(vapply(new_constants, is.numeric, logical(1L)))) {
    stop("all constants must be numeric")
  }
  new_const_attr = c(
    attr(mat_list, "const"),
    new_constants
  )
  if (!good_names(names(new_const_attr))) {
    stop("names of list of constants are malformed")
  }
  attr(mat_list, "const") = new_const_attr
  mat_list
}

define_state.EpiMatrixList = function(mat_list, state_mat_nm_list) {
  state_mat_list = mat_list[state_mat_nm_list]
  if (!all(vapply(state_mat_list, is, logical(1L), "EpiMatrixInput"))) {
    stop("state variables can only be EpiMatrixInput types")
  }

}


#' alter.EpiMatrixListDerived = function(mat_list, ...) {
#'   stop("calls to _alter_ must come before all calls to _derive_")
#' }
#'

#' alter.EpiMatrixList = function(mat_list, ...) {
#'   valid_funcs = getOption("MP_valid_alter_funcs")
#'   new_mats = list(...)
#'   for (i in seq_along(new_mats)) {
#'     mat_i = eval(
#'       new_mats[[i]][[2L]],
#'       c(mat_list, valid_funcs),
#'       enclos = emptyenv()
#'     )
#'     new_mats[[i]] = epi_mat(
#'       names(new_mats[i]),
#'       as.matrix(mat_i)
#'     )
#'   }
#'   structure(
#'     c(mat_list, new_mats),
#'     class = c("EpiMatrixList")
#'   )
#' }

# internal functions ------------------------------------------------

# Construct a struc_expanded Object
#
# @param l list of struc objects
# @param d dimensions of the resulting object
# @return struc_expanded object

struc_expanded = function(l, d) {
  stopifnot(is.recursive(l))
  l = lapply(l, as.struc)
  new('struc_expanded', l = l, dims = d)
}

# Expand and Contract struc and struc_expanded Objects
#
# @param x object to contract or expand
# @rdname expand_contract_struc

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
  if (is.null(dim(x))) {
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

which_wraped_name = function(x) {
  grepl("^[ ]*\\([ ]*[a-zA-Z]{1}[a-zA-Z0-9_]*[ ]*\\)[ ]*$", x)
}

unwrap_paren = function(x) {
  name_candidates = sub(
    "^([ ]*\\(?[ ]*)([a-zA-Z]{1}[a-zA-Z0-9_]*)([ ]*\\)?[ ]*)$",
    "\\2",
    x
  )
  if (!McMasterPandemic:::good_names(unique(name_candidates))) {
    stop("could not extract one name per element")
  }
  name_candidates
}

simple_mult = function(x, y) {
  c(outer(x@v, y@v, paste, sep = " * "))
}

recycle_mult = function(x, y) {
  paste(x@v, y@v, sep = " * ")
}

get_elementwise_dimnames = function(x, y) {
  if (!same_dimnames(x, y)) {
    if (is_1by1(x)) {
      dn = dimnames(y)
    } else if (is_1by1(y)) {
      dn = dimnames(x)
    } else if (is_1byn(x) & is_1byn(y) & same_colnames(x, y)) {
      dn = list("null", colnames(x))
    } else if (is_nby1(x) & is_nby1(y) & same_rownames(x, y)) {
      dn = list(rownames(x), "null")
    } else {
      stop('if dimnames are different, one operand needs to be 1-by-1')
    }
  } else {
    dn = dimnames(x)
  }
  dn
}

elementwise_compatible_dimnames = function(x, y) {
  if (same_dimnames(x, y) | is_1by1(x) | is_1by1(y)) return(TRUE)
  if (is_1byn(x) & is_1byn(y)) return(same_colnames(x, y))
  if (is_nby1(x) & is_nby1(y)) return(same_rownames(x, y))
  return(FALSE)
}

all_elementwise_compatible_dimnames = function(l) {
  for (i in seq_len(length(l) - 1L)) {
    if (!elementwise_compatible_dimnames(l[[i]], l[[i + 1]])) {
      return(FALSE)
    }
  }
  TRUE
}

good_names = function(nms) {
  not_null = !is.null(nms)
  no_dups = isTRUE(!any(duplicated(nms)))
  good_chars = all(grepl("^[a-zA-Z]{1}[a-zA-Z0-9_]*$", nms))
  no_missing = isTRUE(!any(is.na(nms)))
  no_dups & good_chars & no_missing & not_null
}

assert_good_names = function(x) {
  if (!good_names(names(x))) {
    stop("names are not appropriate for this purpose")
  }
}

is_num_or_chr = function(x) is.numeric(x) | is.character(x)

is_epi_scalar = function(x) {
  right_shape = length(x) == 1L
  right_shape & (is_num_or_chr(x))
}

is_epi_vector = function(x) {
  if (is.matrix(x)) {
    right_shape = xor(nrow(x) == 1L, ncol(x) == 1L)
  } else {
    right_shape = !is_epi_scalar(x)
  }
  right_shape & (is_num_or_chr(x))
}

is_epi_row_vector = function(x) {
  (nrow(x) == 1L) & is_epi_vector(x)
}

is_epi_col_vector = function(x) {
  (ncol(x) == 1L) & is_epi_vector(x)
}

is_epi_matrix = function(x) {
  right_shape = is.matrix(x) & !is_epi_vector(x) & !is_epi_scalar(x)
  right_shape & (is_num_or_chr(x))
}

list_if_not_list = function(x) {
  if (!is.list(x)) x = list(x)
  x
}

get_methods_nlist = function(class) {
  generic_name = attr(methods(class = class), "info")$generic
  setNames(lapply(generic_name, get), generic_name)
}

make_elnames = function(x, nm) {
  if (is_epi_scalar(x)) {
    elnames = nm
  } else if (is_epi_row_vector(x)) {
    elnames = nm %_% colnames(x)
  } else if (is_epi_col_vector(x)) {
    elnames = nm %_% rownames(x)
  } else if (is_epi_matrix(x)) {
    elnames = nm %_% expand_names(rownames(x), colnames(x))
  } else {
    stop("input cannot be converted into an epidemiological matrix")
  }
  elnames
}

# internal utility methods ---------------------------------------

#' Functions for Developers
#' @rdname for_dev
#' @keywords internal
#' @export
setMethod("dim", c(x = 'struc_expanded'), function(x) {
  x@dims
})

# internal mathematical methods ----------------------

setMethod('resolve', c(x = 'struc'), function(x) {
  contract_struc(resolve(expand_struc(x)))
})

setMethod('resolve', c(x = 'struc_expanded'), function(x) {
  (x
   %>% slot('l')

   # 1 times anything should just be that thing
   %>% lapply(gsub, pattern = '(\\(1\\)\\s*\\*{1}|\\*{1}\\s*\\(1\\))', replacement = '')

   %>% lapply(trimws)

   # replace blank elements with the 1s that you removed just above
   %>% lapply(function(y) ifelse(y == '', '(1)', y))

   # since we have only products in struc_expanded objects,
   # we can simplify every term with a zero to be exactly zero
   %>% lapply(sub, pattern = ".*\\(0\\).*", replacement = "(0)")

   # get rid of extraneous zeros
   %>% lapply(function(y) {
     y = y[y != "(0)"]
     if (length(y) == 0L) y = "(0)"
     y
   })
   %>% struc_expanded(d = dim(x))
  )
})

#' Functions for Developers
#' @rdname for_dev
#' @keywords internal
#' @export
setMethod("*", c(e1 = 'struc_expanded', e2 = 'struc'), function(e1, e2) {
  l = lapply(e1@l, simple_mult, e2)
  struc_expanded(lapply(l, struc), d = dim(e1))
})

#' Functions for Developers
#' @rdname for_dev
#' @keywords internal
#' @export
setMethod('*', c(e1 = 'struc_expanded', e2 = 'struc_expanded'), function(e1, e2) {
  # TODO: check that dims match
  l = mapply(simple_mult, e1@l, e2@l, SIMPLIFY = FALSE)
  struc_expanded(lapply(l, struc), d = dim(e1))
})

# internal mathematical functions -------------------------------

# to avoid awkward functional programming with binary operators
scal_mult = function(w, x) w * x

struc_to_list_of_scalars = function(x) {
  l = list()
  x = flatten(x)
  for (i in seq_len(size(x))) {
    l[[i]] = x[i]
  }
  l
}

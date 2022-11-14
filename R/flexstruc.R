# classes ------------------------------------------------

#' Class to Represent Scalar, Vector, or Matrix Structure
#'
#' Objects of class \code{\link{struc}} are matrices that can be used in
#' symbolic mathematical expressions for defining various
#' characteristics of compartmental models. Several mathematical operations
#' exist for \code{\link{struc}} objects, including
#' \link[=elementwise]{elementwise sums and products},
#' \link[=matmult]{matrix multiplication},
#' \link[=sumprod]{total sums and products},
#' \link[=rowcolsums]{marginal row and columns sums},
#' \link[=row_col_mult]{multiplication of matrices by rows and columns},
#' \link[=symbtrans]{matrix transpose},
#' and blockwise matrix construction functions (
#'  \code{\link{struc_block}},
#'  \code{\link{struc_stretch}},
#'  \code{\link{struc_bdiag}}
#' ).
#'
#' @md
#' @details
#'
#' Symbolic expressions can be useful for computing complex rates
#' of transition between compartments. For example, the force of
#' infection in models with infectious classes with different severity
#' levels often requires taking a weighted sum of infectious classes,
#' where the weights are stored in a parameter vector.
#'
#' ```{r}
#' severity = c("mild", "severe")
#' I = struc("I" %_% severity)
#' w = struc("w" %_% severity)
#' ```
#'
#' Here the \code{\link{struc}} function was used to construct symbolic
#' math objects, \code{I} and \code{w}, from simple character string
#' vectors. The \code{\link{%_%}} operator \code{\link{paste}}s character
#' objects together with an underscore delimiter.
#'
#' ```{r}
#' "I" %_% severity
#' ````
#'
#' With these \code{I} and \code{w} objects we can programmatically
#' get an expression for the weighted sum.
#'
#' ```{r}
#' sum(w * I)
#' ````
#'
#' This expression could then be used to defines model variables and
#' state transition rates.
#'
#' Once a \code{struc} object is created it can be given
#' \code{dimnames} by ... and creating a \code{struc_dimnames} object.
#' This \code{struc_dimnames} class extends struc by adding names
#' to the rows and columns of the \code{struc} object. For objects
#' with only a single row or single column, it is not necessary to
#' specify names for rows or columns respectively.
#'
#' Continuing the example, we add row names to the column vectors
#' above and take their weighted sum to get the same answer as before.
#'
#' ```{r}
#' rownames(w) = rownames(I) = severity
#' sum(w * I)
#' ````
#'
#' But if we set incompatible row names we get the following error
#'
#' ```{r}
#' rownames(I) = c("wrong", "names")
#' try(sum(w * I))
#' ```
#'
#' In this way, the \code{struc_dimnames} class allows for stronger tests of
#' dimensional compatibility of matrices used in symbolic computations.
#' These tests require that the names of dimensions are also compatible,
#' in addition to the numbers of rows and columns. These stronger
#' tests can be useful for making it more likely that symbolic matrix
#' operations are substantively meaningful, as opposed to just
#' computationally possible. For example, consider a matrix with
#' rows representing different levels of symptom severity, and another
#' with columns representing different levels of vaccination. If the
#' number of severity levels equals the number of vaccination levels,
#' it will be possible to multiply these two matrices. But this
#' multiplication is unlikely to be meaningful, and if the two matrices
#' are stored as \code{struc_dimnames} objects with descriptive dimnames
#' it will correctly not be possible to multiply them together.
#'
#' # Lifecycle
#'
#' The \code{\link{EpiMatrix-class}} is intended to superseed the
#' \code{struc} class. Currently however, the \code{struc} class is
#' the main way to manage model structure.
#'
#' @param e1 first argument of a binary operator
#' @param e2 second argument of a binary operator
#' @param x \code{struc} object as an argument
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
#' @concept symbolic
#' @export
setClass("struc", representation(v = "character", dims = "numeric"))

#' @rdname struc-class
#' @export
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

#' Class to Represent Matrices for use in Compartmental Models
#'
#' This is a new class that is intended to eventually replace
#' \code{\link{struc-class}}. The benefit of this new class is
#' that it is not neccesary to construct special matrices and
#' vectors (see \code{\link{struc-class}}) in order to define matrix
#' algebraic expressions for use in compartmental models. Instead,
#' the user simply supplies plain numerical R vectors and matrices
#' with consistent names and dimnames, and if special symbolic
#' objects are required by the engine they are created automatically
#' without the user needing to be aware of this. As the engine
#' evolves to require less symbolic manipulation of matrices, this
#' new class will work just as well and so it serves to soften the
#' transition to the new engine (\url{https://canmod.net/misc/cpp_side}).
#'
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

#' @rdname EpiMatrix-class
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

#' @rdname EpiMatrix-class
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

# constructors ------------------------------------------------

#' Construct a struc Object
#'
#' Create model structure objects (see \code{\link{struc-class}}), which are
#' symbolic scalars, vectors, and matrices with elements that are expressions
#' involving parameters, and state variables.
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

#' Construct an Epidemiological Matrix
#'
#' Construct a single object of class \code{\link{EpiMatrix-class}}. Typically,
#' the \code{\link{epi_mat_list}} function is used to construct several
#' such epidemiological matrices.
#'
#' @param nm Length-one character vector giving the name of the matrix
#' @param x numeric or character vector or matrix with either names
#' (for vectors) or dimnames (for matrices)
#'
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

#' Epidemiological Matrix List
#'
#' Create a list of matrices to be used as variables in an
#' epidemiological model. This is currently an experimental
#' feature that is planned to superseed \code{\link{struc}}.
#' For more information on how the matrices are stored see
#' \code{\link{EpiMatrix-class}}.
#'
#' @md
#' @details
#'
#' The main benefit of this method of storing lists of matrices,
#' is that the user specifies a list of consistently named
#' numeric vectors and/or matrices that can be manipulated
#' symbolically. For example, here we supply the numeric inputs
#' required to compute the force of infection for a simple model
#' that includes symptomatic status.
#'
#' ```{r}
#' model_vars = epi_mat_list(
#'   beta = c(mild = 0.2, severe = 0.5),
#'   I = c(mild = 50, severe = 2),
#'   N = 1000
#' )
#' ```
#'
#' We can add symbolically derived matrices to this list, without
#' having to create an intermediate \code{\link{struc}} object. We
#' do this by using the \code{\link{derive}} function.
#'
#' ```{r}
#' model_vars = derive(model_vars, foi ~ sum(beta * I) * inverse(N))
#' model_vars$foi
#' ```
#'
#' This symbolic expression for the force of infection is written in
#' terms of the elements of the numeric input matrices and vectors.
#' We can print out these names using the \code{\link[=dims]{names}}
#' function.
#'
#' ```{r}
#' names(model_vars$beta)
#' names(model_vars$I)
#' names(model_vars$N)
#' ```
#'
#' But the \code{epi_mat_list} machinery hides all of these details
#' involved with element names from the user, unless they want to see them.
#' All the user needs to consider are matrix algebraic operations on the
#' objects in the \code{epi_mat_list}.
#'
#' The price for having these details hidden is to construct consistently
#' named objects. If we alter the above example so that the names of
#' the severity classes do not line up, then we will get an error.
#'
#' ```{r}
#' model_vars = epi_mat_list(
#'   beta = c(mild = 0.2, severe = 0.5),
#'   I = c(bad = 50, name = 2),
#'   N = 1000
#' )
#' try(derive(model_vars, foi ~ sum(beta * I) * inverse(N))$foi)
#' ```
#'
#' See \code{\link{epi_mat_names}} for utilities that can help with
#' making consistent naming choices.
#'
#' @param ... named arguments giving numeric or character vectors (with names)
#' or matrices (with dimnames) to be included in an epidemiological model,
#' or named lists of such vectors or matrices. see \code{\link{epi_mat_names}}
#' for utilities that can help with making consistent naming choices.
#'
#' @return An object of class \code{"EpiMatrixList"}, which can be used to
#' simplify the construction of a \code{\link{flexmodel}}
#'
#' @export
epi_mat_list = function(...) {
  mats = unlist(lapply(list(...), list_if_not_list), recursive = FALSE)
  l = mapply(epi_mat, names(mats), mats, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  if (!all(vapply(l, is, logical(1L), class2 = "EpiMatrix"))) {
    return("all matrices must be of class EpiMatrix")
  }
  structure(
    setNames(l, names(mats)),
    class = "EpiMatrixList",
    const = list()
  )
}

epi_const_list = function(...) {
  l = unlist(lapply(list(...), list_if_not_list), recursive = FALSE)
  if (!all(vapply(l, is_num_or_chr, logical(1L)))) {
    return("all constants must be numeric")
  }
  structure(l, class = "EpiConstList")
}

# generic definitions ------------------------------------------------

#' Derive
#'
#' Derive symbolic expressions of existing epidemiological variables, and
#' add them to the list. See \code{\link{epi_mat_list}} for more details
#' on usage.
#'
#' @param mat_list object containing a list of matrix-like objects
#'
#' @export
derive = function(mat_list, ...) {
  UseMethod("derive")
}

const = function(mat_list, ...) {
  UseMethod("const")
}

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

#' Dimnames from Template
#'
#' Make the dimnames of \code{x} be the same as \code{template}
#'
#' @param x object with or without dimnames of the same dimensions as
#' \code{template}
#' @param template object with dimnames to apply to \code{x}
#'
#' @export
setGeneric("dimnames_from_template", function(x, template) {
  standardGeneric("dimnames_from_template")
})

#' Names to Values
#'
#' Replace the values of an object with the names of its elements.
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


#' Names from Values
#'
#' Get element names of an array-like object by processing the
#' values of the object.
#'
#' @param x object from which to create element names from its values
#'
#' @export
names_from_values = function(x) {
  UseMethod("names_from_values")
}

#' Symbolic Complement and Inverse
#'
#' @param x character vector, \code{\link{struc}} object, or
#' \code{\link{EpiMatrixInput-class}} object
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

#' Inner Product
#'
#' Compute the inner product between two objects.
#'
#' @param x first operand to an inner product
#' @param y second operand to an inner product
#' @export
setGeneric("inner", function(x, y, ...) {
  standardGeneric("inner")
})

#' Linear Combination
#'
#' Compute the linear combination of a list of objects
#'
#' @param w object containing weights for the linear combination
#' @param l list of objects to linearly combine, with one object per
#' element in \code{w}
#'
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

#' Extract Values
#'
#' Extract the values of an object with structure.
#'
#' @param x an object from which to extract values
#'
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

#' Symbolic Matrix Subsetting
#'
#' @name subset
#' @concept symbolic
NULL

#' @rdname subset
#' @method `[` struc
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

#' @rdname subset
#' @export
`[.EpiMatrix` = function(x, i, j, ..., drop = FALSE) {
  epi_mat(x@name, x@default[i, j, drop = FALSE])
}

epi_mat_replace_utility = function(x, i, j, value) {
  y = x@default
  y[i, j] = c(value)
  epi_mat(x@name, y)
}

#' @rdname subset
#' @export
`[<-.EpiMatrixInput` = function(x, i, j, ..., value) {
  if (!is.numeric(value)) {
    stop("can only replace values of an epidemiological input matrix with numerical values")
  }
  epi_mat_replace_utility(x, i, j, value)
}

#' @rdname subset
#' @export
`[<-.EpiMatrixDerived` = function(x, i, j, ..., value) {
  if (!is.character(value)) {
    stop("can only replace values of an epidemiological input matrix with numerical values")
  }
  epi_mat_replace_utility(x, i, j, value)
}

#' @rdname subset
#' @export
`[.EpiMatrixList` = function(x, i, j, ..., drop = FALSE) {
  cls = class(x)
  y = unclass(x)[i, drop = FALSE]
  class(y) = cls
  y
}

#' @rdname subset
#' @export
setMethod("diagonal", c(x = "struc"), function(x) {
  struc(diag(as.matrix(x)))
})

#' @export
as.double.EpiMatrixList = function(x, ...) {
  if (length(select_derived_mats(x)) != 0L) {
    warning("list contained derived matrices, which cannot be converted to numeric")
  }
  x = select_input_mats(x)
  setNames(
    unlist(lapply(x, values), use.names = FALSE, recursive = FALSE),
    unlist(lapply(x, names), use.names = FALSE, recursive = FALSE)
  )
}

#' Dimensions for Symbolic Objects
#'
#' @name dims
#' @concept symbolic
NULL

#' @rdname dims
#' @export
setMethod("dim", c(x = 'struc'), function(x) {
  x@dims
})

#' @rdname dims
#' @export
setMethod("dim", c(x = "EpiMatrix"), function(x) dim(x@default))

#' @rdname dims
#' @export
setMethod("dimnames", c(x = "struc_dimnames"), function(x) {
  x@dimnames
})

#' @rdname dims
#' @export
setMethod("names", c(x = "EpiMatrix"), function(x) x@elnames)

#' @rdname dims
#' @export
setMethod("dimnames<-", c(x = "struc", value = "list"), function(x, value) {
  if (is_nby1(x) & is.null(value[[2]])) {
    value[[2]] = ""
  } else if (is_1byn(x) & is.null(value[[1]])) {
    value[[1]] = ""
  }
  new("struc_dimnames", v = x@v, dims = x@dims, dimnames = value)
})

#' @rdname dims
#' @export
setMethod("dimnames", c(x = "EpiMatrix"), function(x) dimnames(x@default))

#' @rdname dims
#' @export
setMethod("dimnames_from_template", c(x = "struc", template = "struc_dimnames"), function(x, template) {
  dimnames(x) = dimnames(template)
  x
})

#' @rdname size
#' @export
setMethod("size", c(x = "struc"), function(x) {
  length(x@v)
})

#' @rdname size
#' @export
setMethod("size", c(x = "EpiMatrix"), function(x) prod(dim(x@default)))

#' @rdname names_to_values
#' @export
names_to_values.numeric = function(x) {
  assert_good_names(x)
  x[] = names(x)
  x
}

#' @rdname names_to_values
#' @export
names_to_values.character = function(x) {
  assert_good_names(x)
  x[] = names(x)
  x
}

#' @rdname names_to_values
#' @export
names_to_values.matrix = function(x) {
  x[] = make_elnames(x, "")
  x
}

#' @rdname names_to_values
#' @export
names_to_values.EpiMatrixDerived = function(x) {
  x[] = names(x)
  x
}

#' @rdname names_to_values
#' @export
names_to_values.EpiMatrixInput = function(x) {
  y = x@default
  y[] = names(x)
  epi_mat(x@name, y)
}

#' @rdname names_to_values
#' @export
names_to_values.array = function(x) stop("this feature is not developed")

#' @rdname names_from_values
#' @export
names_from_values.EpiMatrixDerived = function(x) {
  # could fail, but failure is probably correct as long as the error message
  # is cleaned up. the reason that failure is correct is that not all symbolic
  # matrices have one single name in each cell
  unwrap_paren(as.character(x))
}

#' @rdname names_from_values
#' @export
names_from_values.struc = function(x) {
  # could fail, but that's probably correct as long as the error message
  # is cleaned up
  unwrap_paren(as.character(x))
}

#' @rdname values
#' @export
setMethod("values", c(x = "EpiMatrix"), function(x) unname(flatten(x)))

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

setMethod("show", signature = "struc", function(object) {
  rn = as.character(seq_len(nrow(object)))
  cn = as.character(seq_len(ncol(object)))
  dimnames(object) = list(rn, cn)
  print(object)
})

setMethod("show", signature = "struc_dimnames", function(object) {
  is_one_row = nrow(object) == 1L
  is_one_col = ncol(object) == 1L
  rn = rownames(object)
  cn = colnames(object)
  cat(
    "struc object with ",
    pluralizer(nrow(object), "row"),
    " and ",
    pluralizer(ncol(object), "column"),
    ":\n\n", sep = '')
  for (i in 1:nrow(object)) {
    for (j in 1:ncol(object)) {
      if (!is_one_row) cat('Row', rn[i], '\n')
      if (!is_one_col) cat('Column', cn[j], '\n')
      cat(trimws(gsub('\\+', '\\+\n', as.matrix(object)[i, j])),
          '\n', sep = '')
    }
  }
})

#' @export
setMethod("show", c(object = "EpiMatrix"), function(object) print(object@default))

# exposed utility functions --------------------------

#' Select Epidemiological Matrices
#'
#' Select different types of epidemiological matrices from an
#' \code{\link{EpiMatrixList}} object.
#'
#' @param mat_list an \code{\link{EpiMatrixList}} object.
#' @name select_epi_mats
NULL

#' @describeIn select_epi_mats Select \code{\link{EpiMatrixDerived-class}}
#' objects
#' @export
select_derived_mats = function(mat_list) {
  mat_list[vapply(mat_list, is, logical(1L), class2 = "EpiMatrixDerived")]
}

#' @describeIn select_epi_mats Select \code{\link{EpiMatrixInput-class}}
#' objects
#' @export
select_input_mats = function(mat_list) {
  mat_list[vapply(mat_list, is, logical(1L), class2 = "EpiMatrixInput")]
}

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

#' Add Derived Matrices
#'
#' @param model \code{\link{flexmodel}} object
#' @param mat object of \code{EpiMatrixDerived-class}
#'
#' @return \code{\link{flexmodel}} object with intermediate computations
#' that could be used in calculating rates and/or including in simulation
#' output
#' @export
add_derived_matrix = function(model, mat) {
  vec_factr(model, names(mat), values(mat))
}

#' @param mat_list object of class \code{EpiMatrixList}, containing at least
#' some objects of \code{EpiMatrixDerived-class}
#' @rdname add_derived_matrix
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

#' Helpers for Making Consistent Names and Dimension Names
#'
#' @param x a character vector of names or a numeric or character vector
#' or matrix
#' @name epi_mat_names
NULL

#' @describeIn epi_mat_names create names that include letters,
#' numbers, and underscores only -- this is the recommended format
#' for naming the components of an \code{epi_mat_list}
#' @export
fix_epi_names = function(x) {

  # replace anything except letters, numbers, and underscores
  # with underscores
  y = gsub("[^a-zA-Z0-9_]+", "_", x)

  # remove any repeated underscores
  y = gsub("_{2,}", "_", y)

  # remove trailing and leading underscores
  y = sub("^_{1}", "", y)
  y = sub("_{1}$", "", y)

  # check to see if there are any integers between
  # underscores
  l = strsplit(y, "_")
  is_number = lapply(l, grepl, pattern = "^[0-9]+$")
  uniq_is_number = unique(is_number)

  # if there is a consistent pattern of where the
  # numbers occur in the names, pad them so that
  # all numbers contain the same number of digits
  # -- this is useful for sorting by name
  if (length(uniq_is_number) == 1L) {
    is_number = unlist(uniq_is_number)
    if (any(is_number)) {
      num_digits = max(nchar(unlist(lapply(l, `[`, is_number))))
      fix_numbers = function(v, i) {
        v[i] = formatC(as.integer(v[i]), width = num_digits, flag = "0")
        v
      }
      l = lapply(l, fix_numbers, is_number)
      y = vapply(l, paste0, character(1L), collapse = "_")
    }
  }

  y
}

#' @describeIn epi_mat_names force names to include only digits
#' @export
force_numeric_epi_names = function(x) {
  y = fix_epi_names(gsub("[^0-9_]+", "_", x))
  if (any(nchar(y) == 0L)) {
    y = fix_epi_names(seq_along(x))
  }
  y
}

#' @describeIn epi_mat_names fix the names of an object to be passed
#' to \code{epi_mat_list}
#' @export
set_fixed_epi_names = function(x) setNames(x, fix_epi_names(x))

#' @describeIn epi_mat_names fix the dimension names of an object to be
#' passed to \code{epi_mat_list}
#' @export
set_fixed_epi_dimnames = function(x) {
  if (!is.array(x)) return(set_fixed_epi_names(x))
  dimnames(x) = lapply(dimnames(x), fix_epi_names)
  x
}

#' @describeIn epi_mat_names fix the names to contain only digits of an
#' object to be passed to \code{epi_mat_list}
#' @export
set_numeric_epi_names = function(x) setNames(x, force_numeric_epi_names(x))

#' @describeIn epi_mat_names fix the dimension names to contain only digits of an
#' object to be passed to \code{epi_mat_list}
#' @export
set_numeric_epi_dimnames = function(x) {
  if (!is.array(x)) return(set_numeric_epi_names(x))
  dimnames(x) = lapply(dimnames(x), force_numeric_epi_names)
  x
}

#' @describeIn epi_mat_names sort the elements of a matrix or vector
#' in the order of their dimnames or names
#' @export
sort_epi_dims = function(x) {
  if (is.matrix(x)) {
    return(x[order(rownames(x)), order(colnames(x))])
  } else if (!is.array(x)) {
    return(x[order(names(x))])
  } else {
    stop("invalid input")
  }
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
as.matrix.EpiMatrix = function(x, ...) as(x, "matrix")

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
as.character.EpiMatrix = function(x, ...) as(x, "character")

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

#' Symbolic Elementwise Operators
#'
#' @param e1 Operand as a \code{\link{struc}} object
#' @param e2 Operand as a \code{\link{struc}} object
#'
#' @concept symbolic
#' @name elementwise
NULL

#' @rdname elementwise
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
})

#' @rdname elementwise
#' @export
setMethod("*", c(e1 = "struc_dimnames", e2 = "struc_dimnames"), function(e1, e2) {
  x = as(e1, "struc") * as(e2, "struc")
  dimnames(x) = get_elementwise_dimnames(e1, e2)
  x
})

#' @rdname elementwise
#' @export
setMethod("*", c(e1 = "struc", e2 = "struc_dimnames"), function(e1, e2) {
  elementwise_dimname_harm(e1, e2) * e2
})

#' @rdname elementwise
#' @export
setMethod("*", c(e1 = "struc_dimnames", e2 = "struc"), function(e1, e2) {
  e1 * elementwise_dimname_harm(e2, e1)
})

#' @rdname elementwise
#' @export
setMethod("+", c(e1 = 'struc', e2 = 'struc'), function(e1, e2) {
  if (!same_dims(e1, e2)) {
    stopifnot(is_1by1(e1) | is_1by1(e2))
  }
  pp = paste(e1@v, e2@v, sep = ' + ')
  dim(pp) = dim(e1)
  return(struc(pp))
})

#' @rdname elementwise
#' @export
setMethod("+", c(e1 = "struc_dimnames", e2 = "struc_dimnames"), function(e1, e2) {
  x = as(e1, "struc") + as(e2, "struc")
  dimnames(x) = get_elementwise_dimnames(e1, e2)
  x
})

#' @rdname elementwise
#' @export
setMethod("+", c(e1 = "struc", e2 = "struc_dimnames"), function(e1, e2) {
  elementwise_dimname_harm(e1, e2) + e2
})

#' @rdname elementwise
#' @export
setMethod("+", c(e1 = "struc_dimnames", e2 = "struc"), function(e1, e2) {
  e1 + elementwise_dimname_harm(e2, e1)
})

#' Symbolic Matrix Multiplication
#'
#' @param x Operand as a \code{\link{struc}} object
#' @param y Operand as a \code{\link{struc}} object
#' @param X Operand as a \code{\link{struc}} object
#' @param Y Operand as a \code{\link{struc}} object
#'
#' @concept symbolic
#' @name matmult
NULL

#' @rdname matmult
#' @export
setMethod("%*%", c(x = 'struc', y = 'struc'), function(x, y) {
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
})

#' @rdname matmult
#' @export
setMethod("%*%", c(x = "struc_dimnames", y = "struc_dimnames"), function(x, y) {
  if (!isTRUE(all.equal(colnames(x), rownames(y)))) {
    stop("cannot multiply two matrices with dimnames unless the colnames of the first are the same as the rownames of the second")
  }
  z = as(x, "struc") %*% as(y, "struc")
  dimnames(z) = list(rownames(x), colnames(y))
  z
})

#' @rdname matmult
#' @export
setMethod("%*%", c(x = "struc", y = "struc_dimnames"), function(x, y) {
  if (is_1by1(x)) {
    dimnames(x) = list("", "")
  } else if (!is_1byn(x)) {
    stop("dimnames cannot be automatically harmonized. please be explicit by ensuring that the colnames of x match the rownames of y.")
  }
  colnames(x) = rownames(y)
  x %*% y
})

#' @rdname matmult
#' @export
setMethod("%*%", c(x = "struc_dimnames", y = "struc"), function(x, y) {
  if (is_1by1(y)) {
    dimnames(y) = list("", "")
  } else if (!is_nby1(y)) {
    stop("dimnames cannot be automatically harmonized. please be explicit by ensuring that the colnames of x match the rownames of y.")
  }
  rownames(y) = colnames(x)
  x %*% y
})

#' @rdname matmult
#' @export
setMethod("kronecker", c(X = 'struc', Y = 'struc'), function(X, Y) {
  # TODO: dimname template
  struc_stretch(X, nrow(Y), ncol(Y)) * struc_block(Y, nrow(X), ncol(X))
})

#' Symbolic Matrix Transpose
#'
#' @param x \code{\link{struc}} object
#' @concept symbolic
#' @name symbtrans
#' @export
setMethod("t", c(x = 'struc'), function(x) {
  x@v = c(t(as.matrix(x)))
  x@dims = rev(x@dims)
  x
})

#' @rdname symbtrans
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

#' @rdname flatten
#' @export
setMethod("flatten", c(x = 'struc'), function(x) {
  struc(as.character(x))
})

#' @rdname flatten
#' @export
setMethod("flatten", c(x = "EpiMatrix"), function(x) {
  xx = x@default
  dim(xx) = NULL
  setNames(xx, names(x))
})

#' Summarise Symbolic Expressions as Scalars
#'
#' @concept symbolic
#' @name sumprod
NULL

#' @rdname sumprod
#' @export
setMethod("sum", c(x = 'struc'), function(x, ..., na.rm = FALSE) {
  l = unlist(lapply(c(list(x), list(...)), slot, 'v'))
  struc(paste(l, collapse = ' + '))
})

#' @rdname sumprod
#' @export
setMethod("sum", c(x = "struc_dimnames"), function(x, ..., na.rm = FALSE) {
  z = callNextMethod()
  dimnames(z) = list("", "")
  z
})

#' @rdname sumprod
#' @export
setMethod("prod", c(x = 'struc'), function(x, ..., na.rm = FALSE) {
  # FIXME: does this work right?
  l = unlist(lapply(c(list(x), list(...)), slot, 'v'))
  struc(paste(l, collapse = ' * '))
})

#' @rdname sumprod
#' @export
setMethod("prod", c(x = "struc_dimnames"), function(x, ..., na.rm = FALSE) {
  z = callNextMethod()
  dimnames(z) = list("", "")
  z
})

#' @rdname sumprod
#' @export
setMethod("inner", c(x = 'struc', y = 'struc'), function(x, y) {
  stopifnot(same_dims(x, y))
  sum(x * y)
})

#' Symbolic Row and Column Sums
#'
#' @param x \code{\link{struc}} object
#'
#' @concept symbolic
#' @name rowcolsums
NULL

#' @rdname rowcolsums
#' @export
setMethod("rowSums", c(x = 'struc'), function(x, na.rm = FALSE, dims = 1L) {
  struc(apply(as.matrix(x), 1, function(y) sum(struc(y))@v))
})

#' @rdname rowcolsums
#' @export
setMethod("rowSums", c(x = 'struc_dimnames'), function(x, na.rm = FALSE, dims = 1L) {
  y = callNextMethod()
  dimnames(y) = list(rownames(x), "")
  y
})

#' @rdname rowcolsums
#' @export
setMethod("colSums", c(x = 'struc'), function(x, na.rm = FALSE, dims = 1L) {
  t(struc(apply(as.matrix(x), 2, function(y) sum(struc(y))@v)))
})

#' @rdname rowcolsums
#' @export
setMethod("colSums", c(x = 'struc_dimnames'), function(x, na.rm = FALSE, dims = 1L) {
  y = callNextMethod()
  dimnames(y) = list("", colnames(x))
  y
})

#' Symbolic Linear Combinations
#'
#' @param w Vector of weights as a \code{\link{struc}} object
#' @param l List of \code{\link{struc}} objects of compatible dimensions
#' to be combined linearly
#'
#' @return \code{\link{struc}} object of the same dimensions as those in
#' the list given by the linear combination
#' \code{w[1] * l[[1]] + w[2] * l[[2]] ...}
#'
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

#' @rdname row_col_mult
#' @export
setMethod("row_mult", c(x = "struc", y = "struc"), function(x, y) {
  if (ncol(x) != size(y)) stop("y must be the same size as a row of x")
  z = struc_stretch(flatten(y), nrow(x), 1L) * flatten(x)
  z_values = z@v
  dim(z_values) = dim(x)
  struc(z_values)
})

#' @rdname row_col_mult
#' @export
setMethod("col_mult", c(x = "struc", y = "struc"), function(x, y) {
  if (nrow(x) != size(y)) stop("y must be the same size as a column of x")
  z = struc_block(flatten(y), ncol(x), 1L) * flatten(x)
  z_values = z@v
  dim(z_values) = dim(x)
  struc(z_values)
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

#' @rdname derive
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
      dn = list("", colnames(x))
    } else if (is_nby1(x) & is_nby1(y) & same_rownames(x, y)) {
      dn = list(rownames(x), "")
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

# return version of x that harmonizes cases where y has dimnames, but not x
elementwise_dimname_harm = function(x, y) {
  if (same_dims(x, y)) {
    dimnames(x) = dimnames(y)
  } else if (is_1by1(x)) {
    dimnames(x) = list("", "")
  } else {
    stop("dimnames cannot be harmonized")
  }
  x
}

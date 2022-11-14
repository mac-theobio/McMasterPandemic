#' Fitting Models Sequentially
#'
#' Inputs:
#' 1. observed data splitting function (see \code{\link{make_splitter}})
#' 2. a start date offset
#' 3. time varying parameters
#'
#' @name sequential
NULL

#' Make Splitter
#'
#' Construct a function that will take an observed dataset
#' (see \code{\link{update_observed}}) and return a list of observed
#' datasets. See below for descriptions of the different types of
#' constructors for different types of splitting.
#'
#' @return A function that takes a data frame with a \code{\link{Date}}
#' column called \code{date} and returns the output of a call to the
#' \code{\link{split}} function.
#'
#' @name make_splitter
NULL

#' @describeIn make_splitter Split based on a custom factor vector. The
#' resulting function makes direct use of the \code{\link{split}} function.
#' @param fac Vector that can be coerced to a \code{\link{factor}} to be
#' used to split the \code{data} argument passed to the resulting
#' data splitter.
#' @export
make_splitter_by_fac = function(fac) {
  as.factor(fac)  ## validity checking
  function(data) split(data, fac)
}


#' @describeIn make_splitter Split data into groups with an approximately
#' equal number of rows.
#' @param n_rows Target number of rows in each split of the data.
#' @importFrom lubridate days weeks
#' @export
make_splitter_eq_rows = function(n_rows) {
  force(n_rows)
  function(data) {
    n_rows_data = nrow(data)
    n_splits = ceiling(n_rows_data / n_rows)

    fac = rep(seq_len(n_splits), each = n_rows)

    ## tack on the remainder to the last group of data
    fac = c(fac, rep(max(fac), times = n_rows_data %% n_rows))

    split(data[order(data$date),], fac)
  }
}

#' @describeIn make_splitter Split data into groups with an
#' approximately equal number of days.
#' @param n_days Target number of days in each split of the data.
#' @importFrom lubridate floor_date
#' @export
make_splitter_eq_time = function(n_days) {
  if (packageVersion("lubridate") < package_version("1.9.0.9000")) {
    stop(
      "\nThis data splitter requires lubridate version 1.9.0.9000 or greater.",
      "\nPlease either update your lubridate or use another data splitter",
      "\nconstructor such as make_splitter_eq_rows or make_splitter_by_fac."
    )
  }
  force(n_days)
  function(data) {
    min_date = min(data$date)
    n_days_data = as.integer(diff(range(data$date), units = "days")) + 1L
    grid = min_date + days(seq(from = 0, to = n_days_data, by = n_days))
    fac = floor_date(data$date, grid)
    split(data[order(data$date),], fac)
  }
}

#' @describeIn make_splitter Split data into groups based on custom dates
#' to be used as cut-points
#' @param date_breaks Date vector giving the target first dates in each split
#' of the data.
#' @export
make_splitter_date_breaks = function(date_breaks) {
  if (packageVersion("lubridate") < package_version("1.9.0.9000")) {
    stop(
      "\nThis data splitter requires lubridate version 1.9.0.9000 or greater.",
      "\nPlease either update your lubridate or use another data splitter",
      "\nconstructor such as make_splitter_eq_rows or make_splitter_by_fac."
    )
  }
  date_breaks = as.Date(date_breaks)
  function(data) {
    fac = floor_date(data$date, sort(date_breaks))
    split(data[order(data$date),], fac)
  }
}

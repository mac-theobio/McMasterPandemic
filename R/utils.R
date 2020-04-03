
## attempt convert x to a date unless it already is one
ldmy <- function(x) if (is(x,"Date")) x else lubridate::dmy(x)

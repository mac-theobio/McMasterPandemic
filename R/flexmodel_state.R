#' Get and Set Model State
#'
#' @param model \code{\link{flexmodel}} object
#' @return state vector or data frame of time-variation of parameters
#' @name get_and_set_model_state
NULL

#' @describeIn get_and_set_model_state State vector used to define the model and used as the initial state if \code{do_make_state == FALSE}.
#' @export
state_default = function(model) {
  UseMethod("state_default")
}

#' @describeIn get_and_set_model_state Initial state vector at \code{model$start_date}
#' @export
state_init = function(model) {
  UseMethod("state_init")
}

#' @describeIn get_and_set_model_state Deterministic value of the final state vector at \code{model$end_date}
#' @export
state_final = function(model) {
  UseMethod("state_final")
}

#' @describeIn get_and_set_model_state Random value of the final state vector at \code{model$end_date}
#' @export
state_final_rand = function(model) {
  UseMethod("state_final_rand")
}

#' @describeIn get_and_set_model_state Condensed state vector at \code{model$end_date}.
#' @export
state_final_cond = function(model) {
  UseMethod("state_final_cond")
}

#' @param value named vector of states
#' @rdname get_and_set_model_state
#' @export
`state_default<-` = function(model, value) {
  UseMethod("state_default<-")
}

#' @rdname get_and_set_model_state
#' @export
`state_init<-` = function(model, value) {
  UseMethod("state_init<-")
}


#' @exportS3Method
#' @rdname get_and_set_model_state
state_default.flexmodel = function(model) {
  (model$state
   %>% unlist
   %>% as_vector_no_attr
  )
}

#' @exportS3Method
state_init.flexmodel = function(model) {
  if (!model$do_make_state) return(state_default(model))
  state = (model
    %>% as.flexmodel_one_step()
    %>% simulation_history(add_dates = FALSE)
  )
  (state[1L, names(model$state)]
    %>% unlist
    %>% as_vector_no_attr
  )
}

#' @exportS3Method
state_final.flexmodel = function(model) {
  state = simulation_history(model, add_dates = FALSE)
  (tail(state[names(model$state)], 1L)
    %>% unlist
    %>% as_vector_no_attr
  )
}

#' @exportS3Method
state_final_rand.flexmodel_obs_error = function(model) {
  state = simulation_history(model, add_dates = FALSE, obs_error = TRUE)
  (tail(state[names(model$state)], 1L)
    %>% unlist
    %>% as_vector_no_attr
  )
}

#' @exportS3Method
state_final_cond.flexmodel = function(model) {
  state = simulation_history(model, add_dates = FALSE, condense = TRUE)
  (tail(state, 1L)
    %>% unlist
    %>% as_vector_no_attr
  )
}

#' @export
`state_default<-.flexmodel` = function(model, value) {
  if (!all(names(value) %in% names(model$state))) {
    stop('only currently existing states can be set')
  }
  model$state[names(value)] = value
  model
}

#' @export
`state_init<-.flexmodel` = function(model, value) {
  if (!all(names(value) %in% names(model$state))) {
    stop('only currently existing states can be set')
  }
  model$state[names(value)] = value
  model$do_make_state = FALSE
  model
}

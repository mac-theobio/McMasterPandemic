#' Topological Sort
#'
#' Topologically sort the state names so that states earlier in the flow
#' of individuals come earlier in sorted order.
#'
#' @param model \code{\link{flexmodel}} object
#'
#' @return character vector of state names in topological order
#' @export
topological_sort = function(model) {
  state_order = c()
  rates = model$rates
  state_nms = names(model$state)
  n = length(state_nms)
  for(i in 1:n) {
    if (length(state_nms) == 0L) break
    remaining_inflow = sapply(
      state_nms,
      has_inflow,
      rates,
      state_nms
    )
    state_order = c(state_order, state_nms[!remaining_inflow])
    rates = rates[!unlist(McMasterPandemic:::get_rate_to(rates)) %in% state_order]
    rates = rates[!unlist(McMasterPandemic:::get_rate_from(rates)) %in% state_order]
    state_nms = state_nms[remaining_inflow]
    if (isTRUE(all(remaining_inflow)) & length(state_nms) != 0L & length(remaining_inflow) != 0L) {
      stop("state network is not acyclic (i think), and therefore cannot be topologically sorted")
    }
  }
  state_order
}

has_inflow = function(focal_state, rates, state_nms) {
  stopifnot(focal_state %in% state_nms)
  to_states = (rates
    %>% McMasterPandemic:::get_rate_to()
    %>% unlist
    %>% unique
  )
  focal_state %in% to_states
}

if(FALSE) {
  ## experimenting with
  state_order = base::setdiff(topological_sort(model), c("V", "X"))
  lapply(state_order, has_inflow, model$rates, state_order)
  get_children = function(focal_state, rates) {
    which_parent = unlist(McMasterPandemic:::get_rate_from(rates)) == focal_state
    McMasterPandemic:::get_rate_to(rates)[which_parent]
  }
  get_children("Ip", model$rates)
  e = list2env(list(path_number = 0))
  tree_maker = function(focal_state, e, rates) {
    children = get_children(focal_state, rates)
    if (length(children) == 0L) {
      e$path_number = e$path_number + 1
    }
    path = setNames(lapply(children, tree_maker, e, rates), children)
    if (length(children) == 0L) {
      return(e$path_number)
    } else {
      return(path)
    }
  }
  S_tree = tree_maker("S", e, model$rates)
  S_tree$E$Ia$R
  S_tree$E$Ip$Im$R
  S_tree$E$Ip$Is$H$R
  S_tree$E$Ip$Is$ICUs$H2$R
  S_tree$E$Ip$Is$ICUd$D
  S_tree$E$Ip$Is$D
  S_tree$E$Ip$Is
}

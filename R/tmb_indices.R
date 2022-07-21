# computing indices and data for tmb ----------------------

sum_indices = function(sums, state, params) {
  sum_index_list = lapply(sums, "[[", "sum_indices")
  sumidx = c(seq_len(length(sums))) + length(state) + length(params)
  sumcount = sapply(sum_index_list, length, USE.NAMES = FALSE) %>% unname
  summandidx = c(unlist(sum_index_list, use.names = FALSE))
  clean_if_empty = function(x) {
    if((length(x) == 0L) | is.null(x)) return(integer(0L))
    x
  }
  list(sumidx = clean_if_empty(sumidx),
       sumcount = clean_if_empty(sumcount),
       summandidx = clean_if_empty(summandidx))
}

factr_indices = function(factrs, state_param_sums) {
  if (length(factrs) == 0L) {
    indices = list(
      spi_factr = integer(0L),
      count = integer(0L),
      spi = integer(0L),
      modifier = integer(0L)
    )
    return(indices)
  }
  sp = state_param_sums
  spi_factr = (factrs
    %>% seq_along
    %>% `+`(length(sp))
  )
  spi <- get_var_indx(factrs)
  count <- get_var_counts(factrs)
  modifier <- get_var_modifiers(factrs)
  names(spi) = names(count)  = names(count) = NULL
  indices <- list(
    spi_factr = spi_factr,
    count = count,
    spi = spi,
    modifier = modifier
  )
  return(indices)
}


pow_indices = function(pow, state_param_sum_factr_pow) {
  if (nrow(pow) == 0L) {
    return(init_tmb_indices$pow_indices)
  }
  list(
    powidx = find_vec_indices(pow$pow_nms, state_param_sum_factr_pow),
    powarg1idx = find_vec_indices(pow$pow_arg1_nms, state_param_sum_factr_pow),
    powarg2idx = find_vec_indices(pow$pow_arg2_nms, state_param_sum_factr_pow),
    powconstidx = find_vec_indices(pow$pow_const_nms, state_param_sum_factr_pow)
  )
}


sim_report_expr_indices = function(exprs, init_sim_report_nms) {
  if (length(exprs) == 0L) {
    indices = list(
      sri_output = integer(0L),
      sr_count = integer(0L),
      sri = integer(0L),
      sr_modifier = integer(0L)
    )
    return(indices)
  }
  sri_output = seq_along(exprs) + length(init_sim_report_nms)
  sri <- get_var_indx(exprs)
  sr_count <- get_var_counts(exprs)
  sr_modifier <- get_var_modifiers(exprs)
  names(sri) = names(sr_count)  = names(sr_modifier) = NULL
  indices <- list(
    sri_output = sri_output,
    sr_count = sr_count,
    sri = sri,
    sr_modifier = sr_modifier
  )
  return(indices)
}


lag_diff_indices = function(model) {
  indices_for_one_pattern = function(pattern_input) {
    nms = intermediate_sim_report_names(model)
    indices = grep(pattern_input$var_pattern, nms, perl = TRUE)
    data.frame(
      sri = indices,
      delay_n = rep(pattern_input$delay_n, length(indices))
    )
  }
  (model$lag_diff
    %>% lapply(indices_for_one_pattern)
    %>% Reduce(f = rbind)
    %>% as.list
  )
}


lag_diff_uneven_indices = function(model) {
  n_delay_mat_rows = length(simulation_dates(model, 1L, -1L))
  lag_indices = (model$lag_diff_uneven
    %>% lapply(getElement, 'input_names')
    %>% lapply(find_vec_indices, intermediate_sim_report_names(model))
  )
  if (length(lag_indices) == 0L) {
    return(model)
  }
  lag_breaks = (model$lag_diff_uneven
    %>% lapply(getElement, 'lag_dates')
    %>% lapply(difftime, model$start_date, units = "days")
    %>% lapply(as.integer)
    %>% lapply(`[`, -1L)
    %>% rep(unlist(lapply(lag_indices, length)))
  )
  lag_ns = (model$lag_diff_uneven
    %>% lapply(getElement, 'lag_dates')
    %>% lapply(diff)
    %>% lapply(as.integer)
    %>% rep(unlist(lapply(lag_indices, length)))
  )

  xx = unlist(lag_ns)
  ii = unlist(lag_breaks)
  jj = rep(seq_along(lag_breaks), unlist(lapply(lag_breaks, length)))

  # sparse integer matrices do not seem to be implemented,
  # and the c++ doesn't seem to work without integer matrices,
  # so we are just going to use dense integer matrices
  delay_n = as.matrix(Matrix::sparseMatrix(
    i = ii, j = jj, x = xx,
    dims = c(model$iters + 1L, length(model$lag_diff_uneven))
  ))
  mode(delay_n) = 'integer'
  sri = unlist(lag_indices)
  nlist(sri, delay_n)
}


conv_indices = function(model) {
  indices_for_one_pattern = function(pattern_input) {
    nms = intermediate_sim_report_names(model)
    indices = grep(pattern_input$var_pattern, nms, perl = TRUE)
    conv_par_indices = lapply(pattern_input$conv_pars, find_vec_indices, model$params)
    init_c_delay_cv = model$params[pattern_input$conv_pars$c_delay_cv]
    init_c_delay_mean = model$params[pattern_input$conv_pars$c_delay_mean]
    init_c_prop = model$params[pattern_input$conv_pars$c_prop]
    qmax = length(make_delay_kernel(init_c_prop, init_c_delay_mean, init_c_delay_cv)) + 1L
    #qmax = qgamma(
    #  0.95,
    #  1/init_c_delay_cv^2,
    #  init_c_delay_mean * init_c_delay_cv^2
    #)
    data.frame(
      sri = indices,
      c_prop_idx = rep(conv_par_indices$c_prop, length(indices)),
      c_delay_cv_idx = rep(conv_par_indices$c_delay_cv, length(indices)),
      c_delay_mean_idx = rep(conv_par_indices$c_delay_mean, length(indices)),
      qmax = rep(qmax, length(indices))
    )
  }
  (model$conv
    %>% lapply(indices_for_one_pattern)
    %>% Reduce(f = rbind)
    %>% as.list
  )
}


ratemat_indices <- function(rates, state_params) {
  sp <- state_params
  ratemat_indices <- sapply(rates, `[[`, "ratemat_indices")
  spi <- get_var_indx(rates)
  count <- get_var_counts(rates)
  modifier <- get_var_modifiers(rates)
  names(spi) <- colnames(ratemat_indices) <- names(count) <- NULL
  indices <- list(
    from = ratemat_indices[1, ],
    to = ratemat_indices[2, ],
    count = count,
    spi = spi,
    modifier = modifier
  )
  return(indices)
}

get_var_indx = function(l) {
  {
    lapply(l, function(y) {
      y$factors$var_indx
    }) %>% unlist()
  }
}

get_var_counts = function(l) {
  sapply(l, function(y) {
    nrow(y$factors)
  })
}

get_var_modifiers = function(l) {
  # an 'item' is either a rate or (common) factr
  (l
   %>% unname
   %>% lapply("[[", "factors")
   %>% bind_rows(.id = "item_indx")
   %>% mutate(new_item = as.logical(c(0, diff(as.numeric(item_indx)))))
   %>% mutate(new_prod = as.logical(c(0, diff(as.numeric(prod_indx)))))
   %>% mutate(add = new_prod & (!new_item))
   %>% mutate(modifier = 4 * add + 2 * invrs + compl)
   %>% `$`("modifier")
  )
}

state_mapping_indices = function(
  state,
  eigen_drop_pattern,
  infected_drop_pattern,
  susceptible_pattern) {
  spec_check(introduced_version = "0.1.1",
             feature = "Disease free state updates")

  if(length(susceptible_pattern) == 0L) {
    susceptible_idx = integer(0L)
  } else {
    susceptible_idx = grep(susceptible_pattern, names(state), perl = TRUE)
  }

  eigen = eigen_drop_pattern
  infected = infected_drop_pattern
  all = names(state)
  if((length(eigen) == 0L) | (length(infected) == 0L)) {
    index_set = list(
      x = list(vector('list', 3L)),
      i = list(vector('list', 3L)),
      j = list(vector('list', 3L)))
  } else {
    index_set = make_nested_indices(all, nlist(eigen, infected), invert = TRUE)
  }

  c(
    (index_set
     %>% getElement('i')
     %>% unlist(recursive = FALSE)
     %>% setNames(c("all_to_eigen_idx",
                    "all_to_infected_idx",
                    "eigen_to_infected_idx"))
    ),
    (index_set
     %>% getElement('j')
     %>% unlist(recursive = FALSE)
     %>% setNames(c("all_drop_eigen_idx",
                    "all_drop_infected_idx",
                    "eigen_drop_infected_idx"))
    ),
    nlist(susceptible_idx)
  )
}


disease_free_indices = function(model) {
  unpack(model)

  df_state_par_idx = (disease_free
    %>% lapply(`[`, 'params_to_use')
    %>% unlist(use.names = FALSE)
    %>% find_vec_indices(params)
  )
  df_state_count = (disease_free
    %>% lapply(`[[`, 'states_to_update')
    %>% lapply(length)
    %>% unlist
  )
  df_state_idx = (disease_free
    %>% lapply(`[`, 'states_to_update')
    %>% unlist(use.names = FALSE)
    %>% find_vec_indices(state)
  )
  nlist(df_state_par_idx, df_state_count, df_state_idx)
}

initial_population_indices = function(model) {
  unpack(model)
  infected_idx = which(initial_population$infected == names(params))
  total_idx = which(initial_population$total == names(params))
  nlist(total_idx, infected_idx)
}

initialization_mapping_indices = function(model) {
  unpack(model)
  state_mapping_indices(
    state,
    initialization_mapping$eigen,
    initialization_mapping$infected,
    initialization_mapping$susceptible)
}

linearized_param_indices = function(model) {
  unpack(model)
  lin_param_vals = (linearized_params
                   %>% lapply(`[`, 'update_value')
                   %>% unlist(use.names = FALSE)
  )
  lin_param_count = (linearized_params
                    %>% lapply(`[[`, 'params_to_update')
                    %>% lapply(length)
                    %>% unlist(use.names = FALSE)
  )
  lin_param_idx = (linearized_params
                  %>% lapply(`[`, 'params_to_update')
                  %>% unlist(use.names = FALSE)
                  %>% find_vec_indices(params)
  )

  nlist(lin_param_vals, lin_param_count, lin_param_idx)
}

outflow_indices = function(outflow, ratemat) {
  if(length(outflow) == 0L) {
    return(list(
      row_count = integer(),
      col_count = integer(),
      rows = integer(),
      cols = integer()
    ))
  }
  indices = lapply(outflow, function(o) {
    list(
      state = grep(o$from, rownames(ratemat), perl = TRUE),
      flow = grep(o$to, colnames(ratemat), perl = TRUE)
    )
  })
  n = length(indices)
  indices$row_count = seq(n)
  indices$col_count = seq(n)
  indices$rows = vector()
  indices$cols = vector()

  for(i in 1:n) {
    rows = indices[[i]]$state
    cols = indices[[i]]$flow

    indices$row_count[i] = length(rows)
    indices$col_count[i] = length(cols)

    indices$rows = c(indices$rows, rows)
    indices$cols = c(indices$cols, cols)
  }

  indices

  # TODO: add warning if any state_indices are equal to other state_indices
  #       that are already added to some other call to outflow in the model.
  #       in general state_indices should be mutually exclusive across each
  #       outflow.
  #nlist(state_indices, flow_state_indices)
  #all = names(model$state)
  #lapply(model$outflow, McMasterPandemic::make_nested_indices, x = all)
}

lin_state_timevar_params = function(schedule) {
  stop('function lin_state_timevar_params is under construction')
}


tmb_observed_data = function(model) {

  lp = model$observed$loss_params
  lp_by_var = (lp
   %>% group_by(Variable, Distribution)
   %>% summarise(loss_param_count = length(na.omit(Parameter)))
  )

  lp_for_vars = (lp
    %>% select(-Parameter)
    %>% distinct
    %>% left_join(lp_by_var, by = c("Distribution", "Variable"))
    %>% mutate(loss_id = find_vec_indices(
       Distribution,
       valid_loss_functions
     ))
     %>% mutate(variable_id = find_vec_indices(
       Variable,
       model$condensation_map[final_sim_report_names(model)]
     ))
     %>% select(loss_param_count, loss_id, variable_id)
  )

  lp_for_params = (lp
    %>% filter(!is.na(Parameter))
    %>% mutate(spi_loss_param = find_vec_indices(
       Parameter,
       c(model$state, model$params)
      )
    )
    %>% select(spi_loss_param)
  )

  # initial_table = (model$observed$loss_params
  #  %>% mutate(loss_id = find_vec_indices(
  #    Distribution,
  #    valid_loss_functions
  #  ))
   # %>% mutate(spi_loss_param = find_vec_indices(
   #     Parameter,
   #     c(model$state, model$params, setNames(NA = 0))
   #    )
   #  )
  #  %>% mutate(variable_id = find_vec_indices(
  #    Variable,
  #    model$condensation_map[final_sim_report_names(model)]
  #  ))
  #  %>% select(variable_id, loss_id, spi_loss_param)
  # )
  # variables_by_distributions = (initial_table
  #   %>% group_by(variable_id, loss_id)
  #   %>% summarise(loss_param_count = n())
  #   %>% ungroup
  # )
  # parameters = (initial_table
  #   %>% select(spi_loss_param)
  # )
  comparisons = (model$observed$data
   %>% na.omit
   %>% rename(observed = value)
   %>% mutate(time_step = find_vec_indices(
     as.character(date),
     as.character(simulation_dates(model))
   ))
   %>% mutate(history_col_id = find_vec_indices(
     var,
     model$condensation_map[final_sim_report_names(model)]
   ))
   %>% select(time_step, history_col_id, observed)
   %>% arrange(time_step, history_col_id)
  )

  # HACK: simplify things for now while we only have a
  # single loss function
  # if (nrow(initial_table) == 0L) {
  #   initial_table$loss_param_count = integer(0L)
  # } else {
  #   initial_table$loss_param_count = 1
  # }
  c(
    #as.list(variables_by_distributions),
    #as.list(parameters),
    #as.list(initial_table),
    as.list(lp_for_vars),
    as.list(lp_for_params),
    as.list(comparisons)
  )
}


tmb_opt_params = function(model) {
  indices = init_tmb_indices$opt_params

  if (length(model$opt_params) > 0L) {

    opt_tables = (model$opt_params
      %>% lapply(tmb_opt_form, model$params)
    )
    indices$index_table = (opt_tables
      %>% lapply(getElement, 'd')
      %>% do.call(what = 'rbind')
    )
    indices$hyperparameters = (opt_tables
      %>% lapply(getElement, 'hyperparams_vec')
      %>% unlist
    )
  }
  any_opt_tv_params = sum(unlist(lapply(model$opt_tv_params, length))) > 0L
  if (any_opt_tv_params) {
    opt_tv_tables = index_tv_table = hyperparameters_tv = setNames(
      vector("list", length(model$opt_tv_params)),
      names(model$opt_tv_params)
    )
    # loop over time-varying types (i.e. abs, rel_orig, rel_prev)
    for(tvt in names(opt_tv_tables)) {
      if(length(model$opt_tv_params[[tvt]]) > 0L) {
        opt_tv_tables[[tvt]] = lapply(
          model$opt_tv_params[[tvt]],
          tmb_opt_form,
          model$params,
          model$timevar$piece_wise$schedule,
          tvt
        )
        index_tv_table[[tvt]] = (opt_tv_tables[[tvt]]
          %>% lapply(getElement, 'd')
          %>% do.call(what = 'rbind')
          #%>% arrange(tv_breaks)
        )
        hyperparameters_tv[[tvt]] = (opt_tv_tables[[tvt]]
          %>% lapply(getElement, 'hyperparams_vec')
          %>% unlist
        )
      }
    }
    indices$index_tv_table = bind_rows(index_tv_table)
    indices$hyperparameters_tv = unlist(hyperparameters_tv)
  }
  indices
}

library(ggplot2)
library(McMasterPandemic)
library(tidyr)



make_variant_params<-function(params, number_of_variants){
  param_names = names(params)
  param_names = param_names[param_names!="N"]
  param_names = param_names[param_names!="E0"]
  variant_params = list()
  for(i in 1:number_of_variants){
    new_params = vector(mode = "numeric")
    for(j in param_names){
      new_params=c(new_params, setNames(runif(1), j))
    }
    new_params = c(new_params, N = params[["N"]], E0 = runif(1, min=0, max=params[["N"]]))
    class(new_params) = "params_pansim"
    variant_params = c(variant_params, list(new_params))
  }
  return(variant_params)
}

substater<-function(fixed_status, fixed_strain, variable_status=c("S", "R"), n_strains){
  stopifnot(length(fixed_status)==length(fixed_strain))
  stopifnot(length(fixed_strain)<=n_strains)
  res = McMasterPandemic:::expand_strain_frame(base_states=variable_status, n_strains=n_strains)
  for(i in 1:nrow(res)){
    for(j in 1:length(fixed_strain)){
      res[i, fixed_strain[j]]=fixed_status[j]
    }
  }
  res=unique(res)
  unname(apply(res, 1, paste0, collapse = ""))
}

# make_model<-function(params){
#   model = (flexmodel(
#     params = params,
#     state = make_state(params = params),
#     start_date = "2020-03-10",
#     end_date = "2020-12-10",
#     do_make_state = FALSE,
#   )
#   %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
#   %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
#   %>% add_rate("Ia", "R", ~ (gamma_a))
#   %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
#   %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
#   %>% add_rate("Im", "R", ~ (gamma_m))
#   %>% add_rate("Is", "H", ~
#                  (1 - nonhosp_mort) * (phi1) * (gamma_s))
#   %>% add_rate("Is", "ICUs", ~
#                  (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
#   %>% add_rate("Is", "ICUd", ~
#                  (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
#   %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
#   %>% add_rate("ICUs", "H2", ~ (psi1))
#   %>% add_rate("ICUd", "D", ~ (psi2))
#   %>% add_rate("H2", "R", ~ (psi3))
#   %>% add_rate("H", "R", ~ (rho))
#   %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
#   %>% add_rate("S", "E", ~
#                  (Ia) * (beta0) * (1 / N) * (Ca) +
#                  (Ip) * (beta0) * (1 / N) * (Cp) +
#                  (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
#                  (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
#   %>% add_outflow(".+", "^(S|E|I|H|ICU|D|R)")
#   )
#   return(model)
# }

# make_full_multivariant_model<-function(variant_params, start_date="2020-01-01", end_date="2022-01-01"){
# NoV = length(variant_params) #NoV = number of variants
#
# params = c(N=0)
# #Create a paramaters object for the new model
# for(i in 1:NoV){
#   tmp_params = variant_params[[i]]
#   tmp_names = names(variant_params[[i]])
#   N_loc = which(tmp_names=="N")
#   params["N"] = params["N"]+tmp_params[N_loc]
#   tmp_params = tmp_params[-N_loc]
#   tmp_names = tmp_names[-N_loc]
#   params = c(params, setNames(tmp_params, paste(tmp_names, i, sep="_")))
# }
# class(params) = "params_pansim"
#
# #Create a state object for the new model
# exclusive_state_names = c("S", "E", "Ia", "Ip", "Im", "Is", "R", "H", "H2", "ICUs")
# inclusive_state_names = c("ICUd", "D")
# state_names = exclusive_state_names
# for(i in 1:(NoV-1)){
#   state_names = expand_names(state_names, exclusive_state_names)
# }
# state_names = c(state_names, inclusive_state_names)
# state = rep(0, length(state_names))
# state = setNames(state, state_names)
# for(i in 1:NoV){
#   exposed_subname_vec = c(rep("S", i-1), "E", rep("S", NoV-i))
#   susceptible_subname_vec = rep("S", NoV)
#   exposed_name = paste0(exposed_subname_vec, collapse = "_")
#   susceptible_name=paste0(susceptible_subname_vec, collapse="_")
#   state[exposed_name]=variant_params[[i]]["E0"]
#   if(i==1) state[susceptible_name] = variant_params[[i]]["N"]-variant_params[[i]]["E0"]
#   else state[susceptible_name] = state[susceptible_name]+(variant_params[[i]]["N"]-variant_params[[i]]["E0"])
# }
# model = flexmodel(params=params, state = state, start_date=start_date, end_date=end_date)
#
# return(model)
# }
#
#
# params_orig = read_params("ICU1.csv")
# variant_params = make_full_variant_params(params_orig, 3)
# # variant_models = lapply(variant_params, make_model)
# multivariant_model = make_multivariant_model(variant_params)


make_multivariant_sir_model<-function(params, resistance, recovery, start_date = "2020-01-01", end_date = "2022-01-01"){
  NoV = length(params) # NoV = Number of Variants
  stopifnot(length(resistance)==NoV)
  stopifnot(length(recovery)==NoV)

  params = c(N=0)
  #Create a paramaters object for the new model
  for(i in 1:NoV){
    tmp_params = variant_params[[i]]
    tmp_names = names(variant_params[[i]])
    N_loc = which(tmp_names=="N")
    params["N"] = params["N"]+tmp_params[N_loc]
    tmp_params = tmp_params[-N_loc]
    tmp_names = tmp_names[-N_loc]
    params = c(params, setNames(tmp_params, paste(tmp_names, i, sep="_")))
    params = c(params, setNames(c(resistance[i], recovery[i]), c(paste0("res_", i), paste0("rec_", i))))
  }
  class(params) = "params_pansim"

  state_names = McMasterPandemic:::expand_strain_names(constrained_states = c("I"), n_strains = NoV)


  state = rep(0, length(state_names))
  state = setNames(state, state_names)
  initial_susceptible_state_vec = rep("S", NoV)
  initial_susceptible_state = paste0(initial_susceptible_state_vec, collapse = "")
  state[initial_susceptible_state] = params["N"]
  for(i in 1:NoV){
    initial_infected_state_vec = c(rep("S", i-1), "I", rep("S", NoV-i))
    initial_infected_state = paste0(initial_infected_state_vec, collapse = "")
    state[initial_infected_state] = params[paste0("E0_", i)]
  }
  class(state) = "state_pansim"

model = flexmodel(params = params, state=state, start_date=start_date, end_date=end_date)

# Add force of infection, separately for each strain

for(i in 1:NoV){# i is the strain we are currently adding the force of infection for
  from_states = substater(fixed_status = c("S"), fixed_strain=c(i), n_strains=NoV)
  to_states = substater(fixed_status=c("I"), fixed_strain=c(i), n_strains = NoV)

  beta = vec(rep(paste("beta", i, sep="_"), length(from_states)))
  beta_modifier = vector(mode="character", length=length(from_states))
  for(j in 1:length(from_states)){
      beta_modifier[j] = paste0("res_", (lengths(regmatches(from_states[j], gregexpr("R", from_states[j])))+1))
    }
  beta_modifier=vec(beta_modifier)
  beta_effective = beta * beta_modifier
  invN = vec(rep("1/N", length(from_states)))
  infectious_states = struc(matrix(to_states, ncol=length(to_states), nrow=length(from_states)))
  foi=rowSums((invN%*%t(beta_effective))*infectious_states)
  model = (model%>%vec_rate(from_states, to_states, foi))
  }

return(model)
}






params = c(beta = 0.5, gamma = 0.25, wane_rate = 0.02, N=1000000, E0 = 1000)
variant_params = make_variant_params(params, 3)
multivariant_sir_model = make_multivariant_sir_model(params = variant_params, resistance=c(1, 0.95, 0.80), recovery=c(1, 1.05, 1.15))


exmpl=substater(fixed_status=c("I"), fixed_strain = c(2), n_strain=3)

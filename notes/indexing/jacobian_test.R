library(McMasterPandemic)
library(tools)
library(dplyr)
library(lubridate)

set_spec_version('0.1.1', '../../inst/tmb')

params <- read_params("ICU1.csv")
state <- make_state(params = params)[1:12]
params[['N']] = 1
params[['E0']] = 1e-5
state[] = 0
state[['S']] = 1 - params[['E0']]
state[['E']] =     params[['E0']]

#       S       E      Ia      Ip      Im      Is       H      H2    ICUs    ICUd       D       R
# 0.99999 0.00001 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000

# [[1]]
# [[1]]$state
# [1] "S"
# [[1]]$flow
# [1] "S"
# [[2]]
# [[2]]$state
# [1] "E"    "Ia"   "Ip"   "Im"   "Is"   "H"    "H2"   "ICUs" "ICUd" "D"    "R"
# [[2]]$flow
# [1] "S"    "E"    "Ia"   "Ip"   "Im"   "Is"   "H"    "H2"   "ICUs" "ICUd" "D"    "R"

xx = run_sim_range(params, state, nt = 101,
                   step_args = list(do_hazard = FALSE,
                                    do_exponential = TRUE))

start_date = ymd(20000101)
model <- (init_model(params, state,
                     start_date = start_date,
                     end_date = start_date + days(100),
                     do_hazard = FALSE,
                     do_make_state = FALSE)
          %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
          %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
          %>% add_rate("Ia", "R", ~ (gamma_a))
          %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
          %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
          %>% add_rate("Im", "R", ~ (gamma_m))
          %>% add_rate("Is", "H", ~
                         (1 - nonhosp_mort) * (phi1) * (gamma_s))
          %>% add_rate("Is", "ICUs", ~
                         (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
          %>% add_rate("Is", "ICUd", ~
                         (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
          %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
          %>% add_rate("ICUs", "H2", ~ (psi1))
          %>% add_rate("ICUd", "D", ~ (psi2))
          %>% add_rate("H2", "R", ~ (psi3))
          %>% add_rate("H", "R", ~ (rho))
          # %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
          %>% add_rate("S", "E", ~
                         (Ia) * (beta0) * (1 / N) * (Ca) +
                         (Ip) * (beta0) * (1 / N) * (Cp) +
                         (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
                         (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
          %>% add_outflow("^S$", "^S$")
          %>% add_outflow("^(E|I|H|ICU|D|R)", "^(S|E|I|H|ICU|D|R)")
          %>% update_tmb_indices
)

simulate_changing_ratemat_elements(model)
yy = xx[101,names(state)][2:10]
zz = final_state_vector(model)[2:10]
all.equal(unlist(yy/sum(yy)), c(zz/sum(zz)))



# vector<Type> result(mat.rows());
#
# // Initialize result with zeros
# for (int i=0; i<result.size(); i++)
#   result(i) = 0.0;
result = numeric(12L)

# int startRow = 0;
startRow = 0

# for (int k=0; k<outflow_row_count.size(); k++) { // groups
#   for (int i=startRow; i<startRow+outflow_row_count(k); i++) { // rows in a group
#     int row = outflow_rows[i] - 1;
#
#     int startCol = 0;
#     for (int j=startCol; j<startCol+outflow_col_count(k); j++) { // cols in a row
#       int col = outflow_cols[j] - 1;
#       result[row] += mat.coeff(row, col);
#     }
#     startCol += outflow_cols(k);
#   }
#   startRow += outflow_row_count(k);
# }
for(k in seq_along(model$tmb_indices$outflow$row_count) - 1) {  # groups
  for(i in startRow:(startRow+model$tmb_indices$outflow$row_count[k+1])) { # rows
    print(k %_% i)
  }
}




plot(log(xx$Ia[1:20]))
lines(log(simulate_state_vector(model)$Ia[1:20]))


obj = tmb_fun(model)
obj$fn()
obj$report()
obj$simulate()
final_state_ratio(model)

if(spec_ver_lt('0.1.1')) {
  # no deprecation period for add_parallel_accumulators
  model = add_parallel_accumulators(model, c("X", "N", "P", "V"))
} else {
  model = (model
           %>% add_outflow(".+", "^(S|E|I|H|ICU|D|R)")

           # Update parameters for use with the linearized model
           %>% update_linearized_params('N', 1) # scale population to 1
           %>% update_linearized_params('E0', 1e-5) # perturbation

           # Set the disease-free equilibrium of the linearized model
           %>% update_disease_free_state('S', 'S0') # instead of N

           # Perturb the disease-free equilibrium of the linearized model
           %>% update_disease_free_state('E', 'E0')

           # Define outflows for the linearized model
           %>% add_linearized_outflow("^S$", "^S$")
           %>% add_linearized_outflow("^(E|I|H|ICU|D|R)", "^(S|E|I|H|ICU|D|R)")

           # Define state mappings used to put the initial state values in
           # the correct positions of the initial state vector
           %>% add_state_mappings(

             # regular expression to find states to drop before computing
             # the eigenvector of the linearized system
             eigen_drop_pattern = '^(X|V)',

             # regular expression to find states to drop from the eigenvector
             # before distributing individuals among infected compartments
             infected_drop_pattern = '^(S|D|R)',

             # regular expression to find states in the initial population
             # of susceptibles
             initial_susceptible_pattern = '^S$'
           )

           # Set the total number of individuals and the total number of
           # infected individuals in the initial state vector
           %>% initial_population(total = 'N', infected = 'E0')
  )
}
update_tmb_indices(model)










model <- make_base_model(
  params, state,
  start_date = start_date, end_date = end_date
)
model$do_hazard = FALSE
obj_fun = tmb_fun(model)

obj_fun$fn()
obj_fun$gr()
obj_fun$he()
report <- obj_fun$report()
report$eigenvec
report$j[1, 1]


model$disease_free
model$tmb_indices$disease_free
model$tmb_indices
model$outflow
model$linearized_outflow

model2 = model %>% update_disease_free_state('^I(a|p|m|s)', '^Ca') %>% update_tmb_indices()
names(model2$state)
names(model2$params)
model2$tmb_indices$disease_free$df_state_par_idx
model2$tmb_indices$disease_free$df_state_count
model2$tmb_indices$disease_free$df_state_idx
model2$disease_free$state$simple

model2 = model %>% update_linearized_params('^C(a|p|m|s)', 10) %>% update_tmb_indices()
names(model2$params)
model2$tmb_indices$disease_free$df_param_vals
model2$tmb_indices$disease_free$df_param_count
model2$tmb_indices$disease_free$df_param_idx
model2$disease_free$params



eee = function(m, x, n) {
  for(i in 1:n) {
    x = m %*% x
    x = x / sum(x)
  }
  x
}
c(eee(report$j[2:12, 2:12], model$state[2:12], 102))
yy = eigen(report$j[2:12, 2:12])$vectors[,1]
yy/sum(yy)

zz = c(0.927318,
0.149089,
0.0527015,
0.252728,
0.0107684,
0.00402067,
8.10154e-05,
0.000738624,
0.00059317,
0.000289799,
0.225991)
zz/sum(zz)

xx =
McMasterPandemic:::disease_free_indices(
  model$disease_free,
  model$state, model$params)

McMasterPandemic:::rate_summary(model)
rate_summary(model)

rates_to_keep = (
    (unlist(get_rate_from(model)) %in% model$disease_free$state$eigen$drop_before_eigen)
  | (unlist(get_rate_to(model)) %in% model$disease_free$state$eigen$drop_before_eigen)
) %>% `!` %>% which
c(state[!(names(state) %in% model$disease_free$state$eigen$drop_before_eigen)], params)
tmb_indices(model)
ratemat_indices(model$rates[rates_to_keep], )
model$disease_free$params
model$disease_free$state$simple
model$disease_free$state$eigen

df_idx = disease_free_indices(model$disease_free, c(model$state, model$params))

names(df_idx)

report$j
state

plot(report$V[1,], eigen(report$j)$vectors[1,])

report$V[1,]
eigen(report$j)$vectors[1,]

abline(a = 0, b = 1)
dimnames(jac) = list(new = names(state), old = names(state))
# derivative of the new state (rows) with respect to the old state (cols)
jac
eigen(jac)



rts = rate_summary(model, TRUE)
sratemat = matrix("0", length(state), length(state), dimnames = list(from = names(state), to = names(state)))
sratemat[cbind(rts$from, rts$to)] = rts$formula
flow = McMasterPandemic::struc_block(names(state), col_times = 14, row_times = 1) * struc(sratemat)

struc_jacobian_component = function(vec, new, old) {
  (vec
  %>% resolve
  %>% as.character
  %>% getElement(new)
  %>% struc
  %>% expand_struc
  %>% slot("l")
  %>% getElement(1L)
  %>% as.character
  %>% lapply(strsplit, " \\* ")
  %>% lapply(unlist)
  %>% lapply(function(y) {
    not_dependent = !grepl(paste0('\\(', names(state)[old], '\\)'), y)
    if(all(not_dependent)) y = '(0)'
    y
  })
  %>% lapply(sub, pattern = names(state)[old], replacement = "1")
  %>% lapply(paste0, collapse = " * ")
  %>% unlist
  %>% paste0(collapse = " + ")
  %>% struc
  %>% resolve
  #%>% struc_eval(c(state, params))
  )
}
struc_jacobian_element = function(flow, new, old) {
  c(as.numeric(new == old) +
  struc_jacobian_component(colSums(flow), new, old) -
    struc_jacobian_component(rowSums(flow * p_accum), new, old))
}
struc_jacobian_component(colSums(flow), 1, 1)
struc_jacobian_component(rowSums(flow * p_accum), 1, 1)

p_accum = ('(1)'
  %>% struc_block(14, 12)
  %>% as.matrix
  %>% cbind('(0)')
  %>% cbind('(0)')
  %>% struc
)

non_zero_inds = which(jac != 0, arr.ind = TRUE)

test_jac = function(new, old) {
  all.equal(struc_jacobian_element(flow, new, old), jac[new, old])
}

do.call(mapply, c(list(FUN = test_jac), unname(as.list(expand.grid(1:14, 1:14)))))
expand.grid(1:14, 1:14)[29,]

mm = matrix(0, 2, 4)
for(i in 1:2) {
  for(j in 3:6) {
    mm[i, j - 2] = struc_jacobian_element(flow, i, j)
  }
}


(flow
  %>% rowSums
as.character(expand_struc(struc(as.character(resolve(rowSums(flow)))[new]))@l[[1]]) # outflow
as.character(resolve(colSums(flow)))[new] # inflow

i = 1; j = 1

jac[i, j]
as.character(rowSums(flow))[i] # outflow
as.character(colSums(flow))[i] # inflow
(expand_struc(flow)@l
  %>% lapply(as.character)
  %>% lapply(sub, pattern = ".*\\(0\\).*", replacement = "(0)")
  %>% lapply(struc)
  %>% struc_expanded(c(14, 14))
  %>% contract_struc
  %>% rowSums
)


eigen(jac[1:12, 1:12])

library(iidda)
colSums(flow) %>%
  as.character %>%
  lapply(extract_all_between_paren, contents_pattern = '[^)]*') %>%
  mapply(FUN = `==`, names(state)) %>%
  lapply(any)
flow


 %>% strsplit(' * ') %>% sapply(getElement, 1L)
colSums(flow) %>% as.character %>% strsplit(' * ') %>% sapply(getElement, 1L)
c(struc_eval(rowSums(xx), c(state, params))) - c(struc_eval(colSums(xx), c(state, params)))




r_sim <- run_sim(
  params = params, state = state,
  start_date = start_date,
  end_date = end_date,
  step_args = list(do_hazard = TRUE),
  condense = TRUE,
  use_flex = FALSE
)



# update_disease_free_eigen('^(X|V)', '^(S|D|R)')

nested_chr_indx = function(x, drop_patterns) {
  output = list()
  for(i in seq_along(drop_patterns)) {

  }
}


a = names(state)
index_set = make_nested_indices(a, list(p = '^(X|V)', e = '^(S|D|R)'), invert = TRUE)

with(index_set, x$e)
with(index_set, x$p[i$e$p])


(a_states = names(state))
(i_ap = grep('^(X|V)', a_states, perl = TRUE, invert = TRUE))
(p_states = a_states[i_ap])
# --- not done below
(before_norm_i = grep('^(S|D|R)', nn, perl = TRUE, invert = TRUE))
(oo = nn[jj])
(kk = grep('^(S|D|R)', oo, perl = TRUE, invert = TRUE))
(pp = oo[kk])

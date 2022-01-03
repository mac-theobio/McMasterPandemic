library(McMasterPandemic)
library(dplyr)
library(tidyr)

grep("^MP", names(options()), value = TRUE)
options(MP_flex_spec_version = "0.1.0")
options(macpan_pfun_method = "grep")

start_date <- "2020-02-01"
end_date <- "2020-09-01"

## initialize params
base_params <- read_params("PHAC.csv")
vax_params <- expand_params_vax(
  params = base_params,
  model_type = "twodose"
)
#vax_params[] = lapply(vax_params, `+`, 0.001)
base_state <- make_state(params = base_params)
vax_state <- expand_state_vax(
  x = base_state,
  model_type = "twodose",
  unif = FALSE
)
#vax_state[] = vax_state + 0.001
#vax_sim <- run_sim(
#  params = vax_params,
#  state = vax_state,
#  start_date = start_date,
#  end_date = end_date,
#  condense_args = list(keep_all = TRUE)
#)

M = make_ratemat(vax_state, vax_params, sparse = TRUE)

# organize problem dimensions
(epi_states = c(attr(vax_state, "epi_cat")))
(asymp_cat = c("S", "E", "Ia", "Ip", "R"))
(vax_cat = c(attr(vax_state, "vax_cat")))

`%_%` = function(x, y) paste(x, y, sep = "_")

McMasterPandemic:::pfun_pairs(
  rep(asymp_cat, 2) %_% 'unvax',
  c(asymp_cat, rep("V", 5)) %_% 'vaxdose1',
  M)

rep_rate(
  from = rep(asymp_cat, 2) %_% 'unvax',
  to = c(asymp_cat, rep("V", 5)) %_% 'vaxdose1',
  formula = ~ beta0,
  vax_state, vax_params, NULL, M
)

list(
  data.frame(
    from = asymp_cat %_% 'unvax',
    to   = asymp_cat %_% 'vaxdose1'),
  data.frame(
    from = asymp_cat %_% 'unvax',
    to   = "V"       %_% 'vaxdose1'),
  data.frame(
    from = asymp_cat %_% 'vaxprotect1',
    to   = asymp_cat %_% 'vaxdose2'),
  data.frame(
    from = asymp_cat %_% 'vaxprotect1',
    to   = "V"       %_% 'vaxdose2')
) %>% bind_rows


find_pos_grep()
rate_matrix_lookup(vax_block_dose1) %>% arrange(from_state, to_state)

look = (M
  %>% rate_matrix_lookup
  %>% separate(from_state, c("from_epi_state", "to_vax_cat"), sep = "_")
  %>% separate(to_state, c("to_epi_state", "from_vax_cat"), sep = "_")
  %>% mutate(from_is_asymp = from_epi_state %in% asymp_cat)
  %>% mutate(to_is_asymp = to_epi_state %in% asymp_cat)
)

find_pos = function(reference, target) {
  stopifnot(length(reference) == length(target))
  stopifnot(is.character(reference) & is.character(target))
  setNames(seq_len(length(target)), target)[reference]
}

(look
  %>% filter(from_vax_cat == "unvax",
             to_vax_cat == "vaxdose1")
)

(look
  #%>% filter(is_asymp)
  %>% mutate(is_unvax = vax_cat == 'unvax')
  %>% mutate(is_vaxdose1 = vax_cat == 'vaxdose1')
  %>% filter(is_unvax | is_vaxdose1)
  %>% group_by(epi_state)
  %>% summarise(from_pos = from_pos[is_unvax],
                to_pos = to_pos[is_vaxdose1])
)

# requires that rownames and colnames of the ratemat are identical and in the
# same order
(epi_indices = sort(find_pos_grep(epi_states, )))

add_state_param_sum("asymp_unvax_N",       paste(asymp_cat, "unvax",       sep = "_"))
add_state_param_sum("asymp_vaxprotect1_N", paste(asymp_cat, "vaxprotect1", sep = "_"))

vax_rate = (
  vec(
    "(vax_prop_first_dose) * (1 / asymp_unvax_N)",
    "(1 - vax_prop_first_dose) * (1 / asymp_vaxprotect1_N)"
  )
  * struc("vax_doses_per_day")
)

McMasterPandemic:::pfun_block(asymp_cat, "V", M)
rep_rate(asymp_cat, "V", )


from_regex <- vax_cat[1] ## unvax
to_regex <- vax_cat[2]
ratemat[
  grepl(from_regex, dimnames(ratemat)$from),
  grepl(to_regex, dimnames(ratemat)$to)
] <- vax_block_dose1



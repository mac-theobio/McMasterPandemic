library(ggplot2)
library(McMasterPandemic)
library(tidyr)
library(dplyr)
library(magrittr)


# ### Example "Classic Macpan Ageified Model
# params = read_params("ICU1.csv")
# params_age = expand_params_age(params, beta0 = runif(10))
# state = make_state(params=params_age, type="ICU1")
# state[state!=0]
# sims = (run_sim(params_age, state, condense = FALSE))

params = read_params("ICU1.csv")
beta0 = vector(mode="numeric", length=16)
for(i in 1:16){
  beta0[i]=0.056*i
}

pmat = as.matrix(read.csv("../PremEtAl_contacts_Canada.csv"))
pmat_row_sums = rowSums(pmat)
pmat_normalizer = do.call(cbind, replicate(length(pmat_row_sums), pmat_row_sums, simplify = FALSE ))
pmat = pmat/pmat_normalizer

age_cat = mk_agecats(min=0, max=75, da=5)

dimnames(pmat)=list(age_cat, age_cat)

params_age = expand_params_age(params, age_cat = age_cat, beta0 = beta0, pmat=pmat)

state = make_state(params=params_age)

ageified_model = make_ageified_model(params = params_age, state = state,  start_date="2020-01-01", end_date="2022-01-01", min_age = 0, max_age = 75, do_ageing = TRUE)

ageified_history = simulation_history(ageified_model)

(ageified_history
  %>% select(Date, any_of(names(ageified_model$state)))
  %>% pivot_longer(-Date)
  %>% separate(name, into = c("epi", "age"), sep = "_",remove = FALSE)
  %>% ggplot()
  + facet_grid(epi~age, scales = "free_y")
  + geom_line(aes(Date, value))
)

# plotting_data = data.frame(
#   Date = ageified_history$Date,
#   stringsAsFactors = FALSE
# )
#
# gen_col_names<-function(base_name){
#   name_vec = vector(mode="character", length=16)
#   for(i in 1:16){
#     if(i!=16){
#       lower_age = 5*(i-1)
#       upper_age = 5*i-1
#       age_name = paste(lower_age, upper_age, sep=".")
#     }
#     else age_name = paste0(5*(i-1), ".")
#     name_vec[i]=paste(base_name, age_name, sep="_")
#   }
#   return(name_vec)
# }

# s_vec = gen_col_names("S")
# s_tab = ageified_history[s_vec]
# plotting_data$S = rowSums(s_tab)
#
# e_vec = gen_col_names("E")
# e_tab = ageified_history[e_vec]
# plotting_data$E = rowSums(e_tab)
#
# ia_vec = gen_col_names("Ia")
# ia_tab = ageified_history[ia_vec]
# plotting_data$Ia = rowSums(ia_tab)
#
# ip_vec = gen_col_names("Ip")
# ip_tab = ageified_history[ip_vec]
# plotting_data$Ip = rowSums(ip_tab)
#
# im_vec = gen_col_names("Im")
# im_tab = ageified_history[im_vec]
# plotting_data$Im = rowSums(im_tab)
#
# is_vec = gen_col_names("Is")
# is_tab = ageified_history[is_vec]
# plotting_data$Is = rowSums(is_tab)
#
# h_vec = gen_col_names("H")
# h_tab = ageified_history[h_vec]
# plotting_data$H = rowSums(h_tab)
#
# h2_vec = gen_col_names("H2")
# h2_tab = ageified_history[h2_vec]
# plotting_data$H2 = rowSums(h2_tab)
#
# icus_vec = gen_col_names("ICUs")
# icus_tab = ageified_history[icus_vec]
# plotting_data$ICUs = rowSums(icus_tab)
#
# icud_vec = gen_col_names("ICUd")
# icud_tab = ageified_history[icud_vec]
# plotting_data$ICUd = rowSums(icud_tab)
#
# r_vec = gen_col_names("R")
# r_tab = ageified_history[r_vec]
# plotting_data$R = rowSums(r_tab)
#
# d_vec = gen_col_names("D")
# d_tab = ageified_history[d_vec]
# plotting_data$D = rowSums(d_tab)
#
# (plotting_data
#   %>% pivot_longer(-Date)
#   %>% ggplot()
#   + facet_wrap(~name, scales = "free")
#   + geom_line(aes(Date, value, color=name))
#   )
# (ggplot()+geom_line(data = plotting_data, aes(x=Date, y=S), color="green")
#   +geom_line(data=plotting_data, aes(x=Date, y=E), color="blue")
#   +geom_line(data=plotting_data, aes(x=Date, y=Ia), color="orange1")
#   +geom_line(data=plotting_data, aes(x=Date, y=Ip), color="orange2")
#   +geom_line(data=plotting_data, aes(x=Date, y=Im), color="orange3")
#   +geom_line(data=plotting_data, aes(x=Date, y=Is), color="orange4")
#   +geom_line(data=plotting_data, aes(x=Date, y=H), color="hotpink1")
#   +geom_line(data=plotting_data, aes(x=Date, y=H2), color="hotpink4")
#   +geom_line(data=plotting_data, aes(x=Date, y=ICUs), color="purple1")
#   +geom_line(data=plotting_data, aes(x=Date, y=ICUd), color="purple4")
#   +geom_line(data=plotting_data, aes(x=Date, y=R), color="grey50")
#   +geom_line(data=plotting_data, aes(x=Date, y=D), color="grey20")
#   )
#
# (ggplot()
#   +geom_line(data=plotting_data, aes(x=Date, y=E), color="blue")
#   +geom_line(data=plotting_data, aes(x=Date, y=Ia), color="orange1")
#   +geom_line(data=plotting_data, aes(x=Date, y=Ip), color="orange2")
#   +geom_line(data=plotting_data, aes(x=Date, y=Im), color="orange3")
#   +geom_line(data=plotting_data, aes(x=Date, y=Is), color="orange4")
#   +geom_line(data=plotting_data, aes(x=Date, y=H), color="hotpink1")
#   +geom_line(data=plotting_data, aes(x=Date, y=H2), color="hotpink4")
#   +geom_line(data=plotting_data, aes(x=Date, y=ICUs), color="purple1")
#   +geom_line(data=plotting_data, aes(x=Date, y=ICUd), color="purple4")
#   +geom_line(data=plotting_data, aes(x=Date, y=D), color="grey20")
# )
#
# (ggplot()
#   +geom_line(data=plotting_data, aes(x=Date, y=H), color="hotpink1")
#   +geom_line(data=plotting_data, aes(x=Date, y=H2), color="green")
#   +geom_line(data=plotting_data, aes(x=Date, y=ICUs), color="purple1")
#   +geom_line(data=plotting_data, aes(x=Date, y=ICUd), color="purple4")
#   +geom_line(data=plotting_data, aes(x=Date, y=D), color="grey20")
# )





















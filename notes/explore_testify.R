library(devtools)
load_all("../")

p <- read_params("PHAC_testify.csv")
state <- make_state(params=p)

ratemat <- make_ratemat(state,p)

ratemat_testify <- make_ratemat(state,p,testify=TRUE)

print(ratemat_testify)

beta_vec <- make_betavec(state,p)
beta_vec_testify <- make_betavec(attr(state,"testify"),p)

update_foi(state,p,beta_vec)
update_foi(attr(state,"testify"),p,beta_vec_testify)

## should not work if they are not both the same type 

update_foi(state,p,beta_vec_testify)

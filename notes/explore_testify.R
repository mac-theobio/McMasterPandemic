library(devtools)
load_all("../")

p <- read_params("PHAC_testify.csv")
state<- make_state(params=p)


ratemat <- make_ratemat(state,p)

ratemat_testify <- make_ratemat(state,p,testify=TRUE)

print(ratemat_testify)

## beta_vec <- make_betavec(state,p)
## beta_vec_testify <- make_betavec(attr(state,"testify"),p)

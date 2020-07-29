library(devtools)
load_all("../")

p <- read_params("PHAC_testify.csv")

## Making states and expanding states
state <- make_state(params=p)
state_testified <- expand_stateval(state)


## Making beta_vec wtr states
beta_vec <- make_betavec(state,p)
beta_vec_testified <- make_betavec(state_testified,p)

## Making ratemat
ratemat <- make_ratemat(state,p)
ratemat_testified <- testify(ratemat,p) 

## Updating FOI 
update_foi(state,p,beta_vec)
update_foi(state_testified,p,beta_vec_testified)

## should not work if they are not both the same type 

update_foi(state,p,beta_vec_testified)


### getting run_sim_range working

debug(run_sim)
# debug(run_sim_range)
# debug(do_step)

sim0 <- run_sim(params = p)
sim0_testified_uncondensed <- run_sim(params = p, ratemat_args = list(testify=TRUE),condense = FALSE)

## It is going to fail here where the default is condense = TRUE
sim0_testified_condensed <- run_sim(params = p, ratemat_args = list(testify=TRUE))


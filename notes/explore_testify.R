library(devtools)
load_all("../")

pp <- read_params("PHAC_testify.csv")

## Making states and expanding states
state <- make_state(params=pp)
state_testified <- expand_stateval(state)


## Making beta_vec wtr states
beta_vec <- make_betavec(state,pp)
beta_vec_testified <- make_betavec(state_testified,pp)

## Making ratemat
ratemat <- make_ratemat(state,pp)
ratemat_testified <- testify(ratemat,pp) 

## Updating FOI 
update_foi(state,pp,beta_vec)
update_foi(state_testified,pp,beta_vec_testified)

## should not work if they are not both the same type 

testthat::expect_error(update_foi(state,pp,beta_vec_testified))


### getting run_sim_range working

## debug(run_sim)
# debug(run_sim_range)
# debug(do_step)

sim0 <- run_sim(params = pp)
options(warn=2,error=recover)
## debug(run_sim)
sim0_testified_uncondensed <- run_sim(params = pp, ratemat_args = list(testify=TRUE), condense = FALSE)
## ratemat: 58 x 58
## state: 56
setdiff(colnames(ratemat),names(state))
sim0_testified_condensed <- run_sim(params = pp,
                                    ratemat_args = list(testify=TRUE))

print(head(sim0_testified_condensed))


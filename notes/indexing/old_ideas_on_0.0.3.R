library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)

# run from tests/testthat
spec_version = "0.0.3"
options(MP_flex_spec_version = spec_version)
context(paste0("tests for spec version: ", spec_version))
test_files = "../../inst/tmb/"

params = read_params('ICU1.csv')
state = make_state(params = params)

m = (
  init_model(
    params, state,
    start_date = '2021-08-26', end_date = "2021-09-26",
    timevar_piece_wise = data.frame(
      Date = c("2021-09-01", "2021-09-15", "2021-09-10", "2021-08-28"),
      Symbol = c("mu", "beta0", "beta0", "mu"),
      Value = c(0.5, 0.1, 2, 0.2),
      Type = c("rel_prev", "rel_orig", "rel_prev", "rel_orig")
    ))
  %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
  %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
  %>% add_rate("Ia", "R", ~ (gamma_a))
  %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
  %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
  %>% add_rate("Im", "R", ~ (gamma_m))
  %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
  %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
  %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
  %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
  %>% add_rate("ICUs", "H2", ~ (psi1))
  %>% add_rate("ICUd", "D", ~ (psi2))
  %>% add_rate("H2", "R", ~ (psi3))
  %>% add_rate("H", "R", ~ (rho))
  %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
  %>% add_rate("S",  "E", ~
                 (Ia) * (beta0) * (1/N) * (Ca) +
                 (Ip) * (beta0) * (1/N) * (Cp) +
                 (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
                 (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
  %>% add_parallel_accumulators(c("X", "N", "P", "V"))
  %>% add_tmb_indices()
)

dd = lapply(unname(time_varying_rates(m)), '[[', 'factors') %>%
  bind_rows(.id = 'rate_indx')
ii = with(dd, which(tv))

jj = which(m$tmb_indices$update_ratemat_indices$spi %in% ii)

nbreaks_fac = m$timevar$piece_wise$nbreaks[dd[ii,]$var]
spi_tv_fac = m$timevar$piece_wise$pi_tv_par[dd[ii,]$var] + length(m$state)

iters = 30
spi_tv_fac = c( 2,      7)
nbreaks_fac    = c( 2,      1)
tbreaks    = c(10, 20, 15)
spi = c(20, 21, 50, 24, 52, 1, 53)

b = with(m$timevar$piece_wise$schedule, as.integer(Date - m$start_date))
s = m$timevar$piece_wise$schedule$Symbol
tbreaks = lapply(names(nbreaks_fac), function(v) {b[s == v]}) %>% unlist

ff = function() {

  # index into the spi vector for each time varying factor
  l = c(1L, 3L, 6L, 10L, 14L, 19L)

  # number of time varying factors
  n = length(l)

  # initial time breaks for each time varying factor
  tau = rep(1, n)

  # number of time breaks for each time varying factor
  m = c(3L, 3L, 2L, 2L, 2L, 2L)

  # locations of all time breaks for all time varying factors,
  # organized such that times vary faster than factors and times
  # are all ascending
  TT = c(2L, 6L, 17L, 2L, 6L, 17L, 15L, 20L, 15L, 20L, 15L, 20L, 15L, 20L)

  # indices into the state param vector for all factors, including
  # those that are not time varying
  spi = c(30L, 27L, 30L, 27L, 3L, 15L, 33L, 18L, 4L, 15L, 33L, 19L, 5L,
          15L, 33L, 20L, 36L, 6L, 15L, 33L, 21L, 36L)

  print(spi)
  for(t in 1:30) {
    print(t)
    for(i in 1:n) {
      if((t == TT[tau[i] + sum(m[1:(i-1)])]) & (m[i] >= tau[i])) {
        spi[l[i]] = spi[l[i]] + 1
        tau[i] = tau[i] + 1
        print(tau)
      }
    }
    print(spi)
  }
}

ff = function() {
  spi_tv_fac = ii
  nbreaks_fac = m$timevar$piece_wise$nbreaks[dd[ii,]$var]
  iters = 30
  n_tv_fac = length(spi_tv_fac)
  ivec_tv_fac = rep(1, n_tv_fac)
  spi = m$tmb_indices$update_ratemat_indices$spi
  for(t in 1:iters) {
    print(spi)
    # update `spi` to account for time variation -----
    facs_to_consider = which(ivec_tv_fac <= nbreaks_fac)
    off = 0
    if(length(facs_to_consider) > 0L) {
      for(i_tv_fac in facs_to_consider) {
        #print(i_tv_fac)
        # loop over the time varying factors and
        # check if you are at a break-point. if you are
        # bump the associated spi elements by one
        if(t == tbreaks[ivec_tv_fac[i_tv_fac] + off]) {
          k = spi_tv_fac[i_tv_fac]
          spi[k] = spi[k] + 1
          ivec_tv_fac[i_tv_fac] = ivec_tv_fac[i_tv_fac] + 1
        }
        off = sum(nbreaks_fac[1:(i_tv_fac-1)]) # off + nbreaks_fac[i_tv_fac]
      }
    }

    # update the rate matrix -----
    #...

    # update the state vector -----
    #...
  }
}

sdate <- "2020-02-10"
edate <- "2020-06-01"
time_pars <- data.frame(Date=c("2020-03-10","2020-03-25"),
                        Symbol=c("beta0","beta0"),
                        Value=c(0.5,0.1),
                        Type = c('rel_orig', 'rel_prev'))
restimedep <- run_sim(params,state,start_date=sdate,end_date=edate,
                      params_timevar=time_pars,ndt=20, condense=FALSE)



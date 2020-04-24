## following https://staff.math.su.se/hoehle/blog/2020/04/15/effectiveR0.html
library(R0)
est_rt_wt <- function(ts, GT_obj) {
    end <- length(ts) - (length(GT_obj$GT)-1)
    R0::est.R0.TD(ts, GT=GT_obj, begin=1, end=end, nsim=1000)
}
## what's best GT to use for estimation?  Can we average across 
GT_pmf <- structure( c(0, 0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1), names=0:7)
GT_obj <- R0::generation.time("empirical", val=GT_pmf)
GT_obj

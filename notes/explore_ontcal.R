library(McMasterPandemic)
library(bbmle)
library(dplyr)
library(ggplot2)
library(broom)
library(purrr)
library(dotwhisker)

## utilities
keep_vars <- c("H","death","report")
kv <- . %>% filter(var %in% keep_vars)
gt <- McMasterPandemic:::get_type

## fighting with calibration.
## * Nelder-Mead **appears* to converge, but refitting gives a (much) better answer!
## * profiling and numeric Hessians (numDeriv::hessian) are very badly behaved (exposes lots
##   of rough edges with how profiling that finds new minima are handled
## * decreasing 'reltol' doesn't seem to make much difference

L <- load("../ontario/ontario_calibration_noICU_2brks_prior.RData")
load("../ontario/ontario_clean.RData")
L <- load("alt_cal.RData")
L <- load("alt_cal_tmp.RData")
-logLik(fit1$mle2)  ## original fit
-logLik(fit2$mle2)  ## refitting gives much, much, much better fit
-logLik(fit3$mle2)  ## lower tolerance (doesn't help at all ...?)
-logLik(fit4$mle2)  ## still stuck

## can we make slice1D work to see if these at least *look* like local min? how
##  about numDeriv::hessian() ?

-logLik(pp2)
tt <- (purrr:::map_dfr(list(fit1=fit1$mle2,fit2=fit2$mle2,fit3=pp2),
                tidy,
                .id="fit")
    ## ugh, pp3 gets extra name on params.log_E0
    %>% bind_rows(tibble(term="params.log_E0",estimate=coef(pp2)[["params.log_E0.params_log_E0"]]))
)
pd <- ggstance::position_dodgev(height=0.25)
ggplot(tt,aes(x=estimate,xmin=estimate-2*std.error,
              xmax=estimate+2*std.error,
              y=term,colour=fit)) +
    geom_point(position=pd) +
    geom_linerange(position=pd)

pred1 <- predict(fit1)
pred2 <- predict(fit2)


dd <- ont_all %>% trans_state_vars() %>% kv() %>% gt()
fit <- ont_cal_noICU_2brks_prior
pp <- predict(fit) %>% kv() %>% filter(date>min(dd$date))
gg0 <- ggplot(pp, aes(date,value,colour=var)) + geom_line() + scale_y_log10()
gg0 + geom_point(data=dd) + geom_vline(xintercept=fit$forecast_args$break_dates)

pp2 <- predict(fit,ensemble=TRUE)
plot(pp2)

m <- McMasterPandemic:::mle_fun
p0=coef(fit$mle2)

tt <- tidy(fit$mle2)
dd2 <- fit$mle2@data$data

ggplot(tt,aes(x=estimate,xmin=estimate-2*std.error, xmax=estimate+2*std.error,
              y=term)) + geom_point() + geom_linerange()


f <- fit$forecast_args
f <- f[!names(f) == "fixed_pars"]
tmpf <- function(p=p0) {
    p1 <- p0
    p1[names(p)] <- p
    do.call(m,c(list(p1),fit$mle2@data))
}
tmpf()
vars<-c("params.log_E0","params.log_beta0","logit_rel_beta01","logit_rel_beta02")
tmpf(p0[vars])
-logLik(fit$mle2)


## 4D surface (11025 values)
n <- c(5,5,21,21)
lwr <- rep(0.8,4)
upr <- rep(1.2,4)
x_args <- Map(
    function(x,n,L,U) {
        seq(x*L,x*U, length.out=n)
    },
    p0[vars],n,lwr,upr)

xx <- as_tibble(do.call(expand.grid,x_args))
library(parallel)
options(mc.cores=6)
xxr <- mclapply(seq(nrow(xx)),
                function(i) tmpf(xx[i,]))
xx$res <- unlist(xxr)
xx$sres <- log10(xx$res-min(xx$res,na.rm=TRUE)+0.001)

p1 <- subset(xx,res==min(res,na.rm=TRUE))
gg1 <- (ggplot(xx,
        aes(x=logit_rel_beta01,y=logit_rel_beta02))
    + facet_grid(params.log_beta0 ~ params.log_E0,labeller=label_both)
    + geom_raster(aes(fill=sres))
    + scale_fill_viridis_c()
    + scale_x_continuous(expand=expansion(0,0))
    + scale_y_continuous(expand=expansion(0,0))
    + theme(panel.spacing=grid::unit(0,"lines"))
)
print(gg1
      + geom_point(colour="red",data=p1,size=3)
      + geom_point(colour="blue",data=as_tibble(rbind(p0)),size=3)
      )

xx_center <- dplyr::filter(xx,params.log_E0==median(params.log_E0),
                            params.log_beta0==median(params.log_beta0))
print(gg1 %+% xx_center
      geom_contour(aes(z=sres)))
dplyr::filter(xx_center,logit_rel_beta01==p0[["logit_rel_beta01"]],
              logit_rel_beta02==p0[["logit_rel_beta02"]])


              
## Sobol

## refit starting from 'best-fit' values
new_pars1 <- relist(coef(fit$mle2),f$opt_pars)
new_fit1 <- calibrate(data=dd2, opt_pars=new_pars1, debug_plot=TRUE,
          base_params=f$base_params)

new_fit1$mle2@details  ## 650!!
plot(predict(new_fit1)) +
    geom_line(data=pp2,lty=2) +
    geom_point(data=dd) + geom_vline(xintercept=fit$forecast_args$break_dates)

plot(predict(new_fit1,ensemble=TRUE))

new_pars2 <- relist(coef(fit$mle2),f$opt_pars)
new_fit1 <- calibrate(data=dd2, opt_pars=new_pars1, debug_plot=TRUE,
          base_params=f$base_params)

## can we fit 'from scratch'?


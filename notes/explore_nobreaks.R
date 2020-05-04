## exploring multimodality etc.
library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(viridisLite)
library(tidyverse)
library(anytime)
library(bbmle)

# load(".run_caltest_nobreak.RData")
# load("run_caltest_nobreak.RData")

NLL <- map_dbl(res,~-logLik(.$fit$mle2))
which.max(NLL)
worst <- res[[which.max(NLL)]]

## how bad is it?
print(ggplot(filter(worst$pred,var=="report"),
        aes(date,value))
    + geom_line()
    + scale_y_log10()
    + geom_point(data=worst$simdat)
)

m <- McMasterPandemic:::mle_fun
f <- worst$fit$forecast_args
f <- f[!names(f) == "fixed_pars"]
-logLik(worst$fit$mle2)

dd <- (worst$simdat
    %>% filter(!is.na(value)) 
    %>% mutate(value = round(value))
)


p0=coef(worst$fit$mle2)
tmpf <- function(p=p0) {
    do.call(m,c(list(p,data=dd),f))
}
tmpf()

## constructing a big log-likelihood surface
n <- c(11,11,10)
lwr <- c(2,-0.15,4)
upr <- c(3,-0.05,5)
x_args <- Map(
    function(x,n,L,U) {
        seq(L,U, length.out=n)
    },
    p0,n,lwr,upr)

xx <- as_tibble(do.call(expand.grid,x_args))
library(parallel)
options(mc.cores=1)
xxr <- mclapply(seq(nrow(xx)),
                function(i) tmpf(xx[i,]))
xx$res <- unlist(xxr)
xx$sres <- log10(xx$res-min(xx$res,na.rm=TRUE)+0.001)

p1 <- subset(xx,res==min(res,na.rm=TRUE))
gg1 <- (ggplot(xx,
        aes(x=params.log_beta0,y=params.log_E0))
    + facet_wrap(~log_nb_disp ,labeller=label_both)
    + geom_raster(aes(fill=sres))
    + scale_fill_viridis_c()
    + scale_x_continuous(expand=expansion(0,0))
    + scale_y_continuous(expand=expansion(0,0))
    + theme(panel.spacing=grid::unit(0,"lines"))
)
print(gg1
      + geom_point(colour="red",data=p1,size=3)
      + geom_point(colour="blue",data=as_tibble(rbind(p0)),size=3)
		+ geom_point(colour="black",data=as_tibble(rbind(true_pars)),size=3)
      + ggtitle("full range: blue=MLE fit, red=worst, true=black")
      )
# print(gg1 %+% filter(xx,params.log_E0==median(params.log_E0),
#                      log_nb_disp==median(log_nb_disp))
#       + geom_point(colour="blue",data=as_tibble(rbind(p0)), size=3)
#       + ggtitle("zoom in on middle panel")
# )
# 
# n2 <- c(1,31,31,1)
# lwr2 <- c(1,0.97,0.97,1)
# upr2 <- c(1,1.03,1.03,1)
# x2_args <- Map(
#     function(x,n,L,U) {
#         seq(x*L,x*U, length.out=n)
#     },
#     p0,n2,lwr2,upr2)
# xx2 <- as_tibble(do.call(expand.grid,x2_args))
# xxr2 <- mclapply(seq(nrow(xx2)),
#                 function(i) tmpf(xx2[i,]))
# xx2$res <- unlist(xxr2)
# xx2$sres <- log10(xx2$res-max(xx2$res,na.rm=TRUE)+0.001)
# print(gg1 %+% xx2 + geom_contour(aes(z=sres))
#       + geom_point(colour="blue",data=as_tibble(rbind(p0)),size=3)
#       + ggtitle(sprintf("zoom in further on MLE fit ('bad', NLL=%1.2f)",
#                         -round(c(logLik(worst$fit$mle2),2))))
# )
# 
# x3_args <- Map(
#     function(x,n,L,U) {
#         seq(x*L,x*U, length.out=n)
#     },
#     p1[,seq_along(p0)],n2,lwr2,upr2)
# xx3 <- as_tibble(do.call(expand.grid,x3_args))
# xxr3 <- mclapply(seq(nrow(xx3)),
#                 function(i) tmpf(xx3[i,]))
# xx3$res <- unlist(xxr3)
# xx3$sres <- log10(xx3$res-max(xx3$res,na.rm=TRUE)+0.001)
# print(gg1 %+% xx3 + geom_contour(aes(z=sres))
#       + geom_point(colour="red",data=p1,size=3)
#       + ggtitle(sprintf("zoom in further on worst fit ('better', NLL=%1.2f)",
#                         round(max(xx3$res,na.rm=TRUE),2)))
#       )

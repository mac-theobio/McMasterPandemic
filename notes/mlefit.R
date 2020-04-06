library(McMasterPandemic)
library(bbmle)
library(plyr)
library(dplyr) ## second!
library(ggplot2); theme_set(theme_bw())

p <- read_params(system.file("params","ICU1.csv",
                             package="McMasterPandemic"))
    
## n.b. extra factor levels cause problems below
## FIXME: drop levels upstream?
clean_data <- dplyr::as_tibble(
                         droplevels(readRDS("../data/01042020_clean.rds")))

## fit only after March 18 (shortcut to estimating effects of control measures)
fit_data <- dplyr::filter(clean_data,date>=as.Date("2020-03-18"))

do_profile <- TRUE

## assumes equal dispersion for all three data types
## (H, ICU, deaths); is this plausible?
t_fit1 <- system.time(m1 <- mle2(value~dnbinom(
               mu=predfun(exp(log_beta0),exp(log_E0),
                          data=clean_data, params=pp),
               size=exp(logk)),
               data=clean_data,
               start=list(log_beta0=0,log_E0=log(300),logk=1),
               method="Nelder-Mead")
               )
print(t_fit1)

## ~10 seconds (158 function evaluations)

## would like to generalize.  How can I pass in new parameters?? Global??
## mle2 is doing ugly things with evaluation environment ...
fitfun <- function(dist="nbinom",
                   parameters=NULL,
                   start=list(log_beta0=0,log_E0=log(300),logk=1),
                   method="Nelder-Mead",
                   xdata=fit_data,
                   control=list(),
                   params=pp
                   ) {
    form <- switch(dist,
                   nbinom=value~dnbinom(mu=predfun(exp(log_beta0),exp(log_E0),
                                                   data=data, params=params),
                                        size=exp(logk)),
                   stop("unknown distribution")
                   )
    tt <- system.time(
        m <- mle2(form,parameters=parameters,
                  start=start,method=method,control=control,
                  data=data))
    attr(m,"time") <- tt
    return(m)
}

## m1 <- fitfun()

## would like to use update() but doesn't work
## only change is parameters=list(logk ~ var -1 ),
##  i.e. variable specific dispersion params
## FIXME: data is entered in two different places

t_fit1B <- system.time(m1B <- mle2(value~dnbinom(
               mu=predfun(exp(log_beta0),exp(log_E0),
                          data=fit_data, params=pp),
               size=exp(logk)),
               data=fit_data,
               parameters=list(logk ~ var -1 ),
               start=list(log_beta0=0,log_E0=log(300),logk=1),
               method="Nelder-Mead",
               control=list(maxit=2000)))
print(t_fit1B)  ## 40 seconds

t_fit1B <- system.time(m1B <- mle2(value~dnbinom(
               mu=predfun(exp(log_beta0),exp(log_E0),
                          data=fit_data, params=pp),
               size=exp(logk)),
               data=fit_data,
               parameters=list(logk ~ var -1 ),
               start=list(log_beta0=0,log_E0=log(300),logk=1),
               method="Nelder-Mead",
               control=list(maxit=2000)))


t_fit2 <- system.time(m2 <- mle2(value~dlnorm(
                          meanlog=predfun(exp(log_beta0),exp(log_E0),
                                          data=clean_data, params=pp),
                          sdlog=exp(log_sdlog)),
                          data=clean_data,
               ## var-specific log std devs
               parameters=list(log_sdlog~var-1),
               start=list(log_beta0=0,log_E0=log(300),
                          log_sdlog=0),
               method="Nelder-Mead")
               )
print(t_fit2) ## 10 seconds

sqrt(diag(vcov(m1))) ## sd, on log scale (proportional sd)
nb0prof <- NULL
if (do_profile) {
    t_prof1 <- system.time(nb0prof <- profile(m1,trace=TRUE))
    ## 99 seconds
    confint(nb0prof)
}
## Wald and profile CIs are v. similar (for log-scaled params)
confint(m1,method="quad")
exp(confint(m1,method="quad"))
cc <- coef(m1)

## predict for specific coefs
pf <- function(cc,...) {
    predfun(exp(cc[["log_beta0"]]),exp(cc[["log_E0"]]),
            clean_data,
            pp,
            ...)
}
pred_full <- pf(coef(m1), values_only=FALSE)
##

## FIXME: allow for longer time vector
get_preds <- function(fit,seed=NULL,level=0.8) {
    if (!is.null(seed)) set.seed(seed)
    pars <- MASS::mvrnorm(n=200,mu=coef(fit),Sigma=vcov(fit))
    all_pred <- plyr::adply(pars,
                            .margins=1, ## rowwise
                            .fun=pf,
                            values_only=FALSE,
                            .progress="text")
    ## ~ 10 seconds
    preds <- (all_pred
        %>% as_tibble()
        %>% dplyr::select(-X1)
        %>% group_by(date,var)
        %>% dplyr::summarise(median=quantile(value,0.5),
                      lwr=quantile(value,(1-level)/2),
                      upr=quantile(value,(1+level)/2))
        %>% full_join(pf(coef(fit),values_only=FALSE),by=c("date","var"))
    )
    return(preds)
}

fitList <- list(nb0=m1,nb_vardisp=m1B,ln_vardisp=m2)
## m2 has non-pos-def vcov? ??
predList <- purrr::map(fitList[1:2], get_preds)

save(fitList,predList,file="mlefit.RData")
## straight NB fit is terrible
gg0 <- (ggplot(predList[["nb_vardisp"]],aes(date,colour=var,fill=var))
    + geom_line(aes(y=value))
    + geom_ribbon(colour=NA,alpha=0.2,aes(ymin=lwr,ymax=upr))
    + geom_line(aes(y=median),lty=2)
    + geom_point(data=clean_data,aes(y=value))
)
gg0 + scale_y_log10()

cc <- coef(m1B)
pp2 <- pp
pp2[["beta0"]] <- exp(cc[["log_beta0"]])
pp2[["E0"]] <- exp(cc[["log_E0"]])
write_params(pp2,"PARAMFILE.csv","NBdispvar fit to H/D/ICU >= March 18")
## thrown off by noise in the first few values?
## control measures/nonlinearity?

library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr); theme_set(theme_bw())

## if we want to do more spliney stuff the splines2 package seems useful:
## https://cran.r-project.org/web/packages/splines2/vignettes/splines2-intro.html

source("makestuff/makeRfuns.R")

makeGraphics()

tmax <- 140

fitmax <- 126 # (note fitmax<=tmax)

t <- 1:tmax 
r <- 0.04
h <- 2
R0 <- 3
c <- 2

Q <- exp(r*t)
S <- 1/(1+Q)
I <- Q/(1+Q^2)
logRt <- log(R0*S^h*I^(-c))

Rf <- data.frame(
	t=t
	, logRt=ifelse(t>fitmax, NA, logRt)
)
print(Rf)

## ???
## default na.action (options("na.action")) is na.omit(), so I have
##  no idea why these should give different results ...
form <- logRt~bs(t,df=7)
coef(lm(form,data=na.omit(Rf)))
coef(lm(form,data=Rf))

## model frames (the stage at which NAs are supposed to be discarded)
## are the same size ...
dim(mf1 <- model.frame(form,data=Rf))
dim(mf2 <- model.frame(form,data=na.omit(Rf)))

## but the computed bases are quite different (!)
par(mfrow=c(1,2), las=1)
matplot(as.matrix(mf1[,-1]),ylim=c(0,1),type="l",lty=1,
        main="full data", col=palette())
matplot(as.matrix(mf2[,-1]),ylim=c(0,1),type="l",lty=1,
        main="na.omit()", col=palette())

## conclusion: I don't know why this happens (my naive assumption
## would be that R should Do The Right Thing and give the same answer
## in either case), but it seems safer to explicitly na.omit()
mod_bs <- lm(logRt~bs(t,df=7),data=na.omit(Rf))
co_bs <- coef(mod_bs) 
print(co_bs)

mod_ns <- lm(logRt~ns(t,df=7),data=na.omit(Rf))
co_ns <- coef(mod_ns) 
print(co_ns)

Rf <- (Rf
    %>% mutate(Rt_predict_bs = predict(mod_bs, newdata=Rf),
               Rt_predict_ns = predict(mod_ns, newdata=Rf))
)
## now we get warning messages about ill-conditioned
## bases (but only from the the bs() prediction)

## to construct the logRt from lm:

## this object has the information needed to reconstruct the correct basis
print(pv <- attr(terms(mod_ns),"predvar"))
## it is a language object, so to use it we need to evaluate it, either
eval(pv)[[2]]  ## pv evaluates to a list, we take the second element
## or
eval(pv[[3]])  ## ns(...) is the *third* element of pv (the first is 'list()')

bns <- eval(pv[[3]])
## bns <- ns(t,df=7)
bns <- cbind(1, bns)  ## cbind() automatically converts to a matrix

## we don't need to explicitly convert coef vector to a matrix
## (R is sloppy but effective)
aliP_ns <- bns %*% co_ns

## can do the same for the bs fit (left as an exercise): this version
## is probably wrong
b <- bs(t,df=7)
b <- as.matrix(cbind(rep(1, nrow(b)), b))
aliP_bs <- b %*% matrix(co_bs,ncol=1)

print(ggplot(Rf)
	+ aes(t, logRt)
	+ geom_line()
	+ geom_point(aes(t,Rt_predict_bs), col="red", shape=1)
	+ geom_point(aes(t, aliP_bs), col="blue", shape=3)
	+ geom_point(aes(t, aliP_ns), col="green", shape=5)
)


Rf_long <- tidyr::pivot_longer(Rf,-t)
print(ggplot(Rf_long,aes(t,value, colour=name)) 
      + geom_line(alpha=0.5)
      + geom_point(aes(shape=name),alpha=0.5)
      + scale_shape_manual(values=c(1,3,5))
      + geom_vline(xintercept=fitmax, lty=2)
)

Rf_longdiff <- (tidyr::pivot_longer(Rf,-c(t,logRt))
    %>% mutate(value=value-logRt)
)

print(ggplot(Rf_longdiff,aes(t,value, colour=name)) 
      + geom_line(alpha=0.5)
      + geom_point(aes(shape=name),alpha=0.5)
      + scale_shape_manual(values=c(1,3,5))
      + geom_vline(xintercept=fitmax, lty=2)
      + labs(y="difference from logRt")
)


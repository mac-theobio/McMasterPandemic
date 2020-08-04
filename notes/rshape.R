library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr); theme_set(theme_bw())

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

mod_bs <- lm(logRt~bs(t,df=7),data=Rf)
co_bs <- coef(mod_bs) 
print(co_bs)

mod_ns <- lm(logRt~ns(t,df=7),data=Rf)
co_ns <- coef(mod_ns) 
print(co_ns)

Rf <- (Rf
	%>% mutate(Rt_predict_bs = predict(mod_bs, newdata=Rf))
)

#to construct the logRt from lm:
b <- bs(t,df=7)
b <- as.matrix(cbind(rep(1, nrow(b)), b))
aliP_bs <- b %*% matrix(co_bs,ncol=1)

bns <- ns(t,df=7)
bns <- as.matrix(cbind(rep(1, nrow(bns)), bns)) #insert col of 1 for intercept
aliP_ns <- bns %*% matrix(co_ns,ncol=1)

print(ggplot(Rf)
	+ aes(t, logRt)
	+ geom_line()
	+ geom_point(aes(t,Rt_predict_bs), col="red", shape=1)
	+ geom_point(aes(t, aliP_bs), col="blue", shape=3)
	+ geom_point(aes(t, aliP_ns), col="green", shape=5)
)


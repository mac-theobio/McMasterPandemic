library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr); theme_set(theme_bw())

source("makestuff/makeRfuns.R")

makeGraphics()

tmax <- 112
fitmax <- 98

t <- 1:tmax #it was 140
r <- 0.04
h <- 2
R0 <- 3
c <- 2

Q <- exp(r*t)
S <- 1/(1+Q)
I <- Q/(1+Q^2)
Rt <- R0*S^h*I^(-c)

Rf <- data.frame(
	t=t
	, Rt=ifelse(t>fitmax, NA, Rt)
)

print(Rf)

mod <- lm(Rt~bs(t,df=7),data=Rf)
co <- coef(mod) 
print(co)

Rf <- (Rf
	%>% mutate(Rt_predict = predict(mod, newdata=Rf))
)

#to construct the Rt from lm:
b <- bs(t,df=7)
b <- as.matrix(cbind(rep(1, nrow(b)), b[,]))
aliP <- b %*% matrix(co,ncol=1)

print(ggplot(Rf)
	+ aes(t, Rt)
	+ geom_line()
	+ geom_point(aes(t,Rt_predict), col="red", shape=1)
	+ geom_point(aes(t, aliP), col="blue", shape=3)
)

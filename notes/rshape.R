library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr); theme_set(theme_bw())

source("makestuff/makeRfuns.R")

fitmax <- 98
t <- 1:122 #it was 140
r <- 0.04
h <- 2
R0 <- 3
c <- 2

Q <- exp(r*t)
S <- 1/(1+Q)
I <- Q/(1+Q^2)

Rt <- R0*S^h*I^(-c)

Rf <- data.frame(
	t=t, Rt=Rt
)

Rf <- (Rf
	%>% mutate (Rt=ifelse(t>fitmax, NA, Rt))
)

print(Rf)

pts <- (ggplot(Rf)
	+ aes(t, Rt)
	+ geom_point()
	+ xlim(range(t))
	+ ylim(c(0, NA))
)

print(pts)

print(pts 
	+ geom_smooth(method=lm, formula = y~ns(x, 7))
	+ ggtitle("natural spline")
)

print(pts 
	+ geom_smooth(method=lm, formula = y~bs(x, 7))
	+ ggtitle("B spline")
)

mod <- lm(Rt~bs(t,7),data=Rf)
mod_nointercept <- update(mod, .~.-1)  
print(coef(mod_nointercept) )


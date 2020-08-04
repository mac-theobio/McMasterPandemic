library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr); theme_set(theme_bw())

source("makestuff/makeRfuns.R")

makeGraphics()

fitmax <- 98 #98
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





mod <- lm(Rt~bs(t,df=7),data=Rf)
cmod <- coef(mod) 
print(cmod)

t_sub <- t[t %in% c(1:fitmax)]

Rf <- (Rf
       %>% mutate (Rt_predict=ifelse(t>fitmax, NA, predict(mod, data=t_sub ) ))
)

pts <- (ggplot(Rf)
        + aes(t, Rt)
        + geom_point()
        +geom_line(aes(t,Rt_predict), col="red")
        + xlim(range(t))
        + ylim(c(0, NA))
)
print(pts) 
 

#to construct the Rt from lm:
b <- bs(t,df=7)
length(b[,1])
plot(t, as.matrix(cbind(rep(1,length(t)) ,b[,])) %*% matrix(cmod,ncol=1))

# no intercept
strmod_nointercept <- update(mod, .~.-1)  
print(coef(mod_nointercept) )



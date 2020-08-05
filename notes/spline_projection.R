# Run spline_fit first
library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

tmax <- 200

t <- 1:tmax 

newdat <- data.frame(t=t)

# Step 1: taking the basis functions for shorter t=fitmax in spline_fit.R and extend the priction
pframe <- (newdat
	%>% mutate(bs_projection = predict(mod_bs,newdata=newdat)
		, ns_projection = predict(mod_ns,newdata=newdat)	)
	%>% gather(key="spline_type",value="pred",-t)
)

gg <- (ggplot(pframe, aes(x=t,y=pred,color=spline_type))
	+ geom_line()
	+ geom_vline(aes(xintercept = 126))
	+ scale_color_manual(values=c("red","blue"))
)

print(gg)

# Step 2: What does it really look like?

t <- 1:tmax
r <- 0.04
h <- 2
R0 <- 3
c <- 2

Q <- exp(r*t)
S <- 1/(1+Q)
I <- Q/(1+Q^2)
logRt <- log(R0*S^h*I^(-c))

Rfnew <- data.frame(t, logRt)

print(gg2 <- gg + geom_point(data=Rfnew, aes(x=t, y=logRt), color="black", alpha=0.3))

## Step 3: Freezing

fitmax <- 126
tmax <- 200
t <- 1:fitmax
#coefficients of bs on shorter time fitmax
co_bs <- coef(mod_bs) 

b <- bs(t,df=7)
b <- as.matrix(cbind(1, b))
bfreeze <- rbind(b,b[rep(126,tmax-fitmax),])

matplot(as.matrix(bfreeze),ylim=c(0,1),type="l",lty=1,
        main="full data", col=palette())

bsfreeze <- bfreeze %*% matrix(co_bs,ncol=1)

t <- 1:tmax

co_bs <- coef(mod_bs)

b <- bs(t,df=7)
b <- as.matrix(cbind(1, b)) 
bscont <- b %*% matrix(co_bs,ncol=1)

aliframe <- data.frame(bsfreeze,bscont,t=1:tmax)

print(gg2 + geom_point(data=aliframe, aes(x=t,y=bsfreeze), color="red", alpha=0.3))
print(gg2 + geom_point(data=aliframe, aes(x=t,y=bscont), color="red", alpha=0.3))


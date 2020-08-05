# Run spline_fit first
library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

tmax <- 140
fitmax <- 126
t <- 1:tmax 

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
  , logRt_all = logRt
)

# Step 1:
# model fits

mod_bs <- lm(logRt~bs(t,df=7),data=na.omit(Rf))
co_bs <- coef(mod_bs) 

mod_ns <- lm(logRt~ns(t,df=7),data=na.omit(Rf))
co_ns <- coef(mod_ns) 

mod_bs <- model.frame(mod_bs)
mod_ns <- model.frame(mod_ns)

print(pv_ns <- attr(terms(mod_ns),"predvar"))
bns <- eval(pv_ns[[3]])
bns <- cbind(1, bns)
# projection for ns
proj_ns <- bns %*% co_ns

print(pv_bs <- attr(terms(mod_bs),"predvar"))
bbs <- eval(pv_bs[[3]])
bbs <- cbind(1, bbs)
proj_bs <- bbs %*% co_bs

Rf <- (Rf 
  %>%
  mutate(ns_projection = proj_ns,
         bs_projection = proj_bs)
  %>% gather(key="spline_type",value="pred",-t,-logRt,-logRt_all)
)


gg <- (ggplot(Rf, aes(x=t,y=pred,color=spline_type))
       + geom_line()
       + geom_vline(aes(xintercept = fitmax))
       + scale_color_manual(values=c("red","blue"))
       + geom_point(aes(x=t,y=logRt_all),col="black",alpha=0.3)
)

print(gg)

Rf <- ( Rf
        %>%
          mutate(difference= pred - logRt_all)
)

print(
  ggplot(Rf, aes(x=t,y=difference,color=spline_type))
  + geom_point(aes(shape=spline_type),alpha=0.5)
  + geom_vline(aes(xintercept = fitmax))
  + scale_color_manual(values=c("red","blue"))
  + labs(y="difference from logRt")
)





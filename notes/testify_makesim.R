## step by step ...

library(McMasterPandemic)
pp <- read_params("PHAC_testify.csv")

## rExp "run an exponential simulation"
s1 <- rExp(pp, testify=TRUE, return_val="sim")
s2 <- rExp(pp, return_val="sim")
stopifnot(identical(s1,s2))  ## has_testing() works!
x <- s1[,-1] ## leave out time
is_S <- grepl("^S",names(x))
matplot(s1$t,x,log="y",type="l",col=1+as.numeric(is_S),
        lty=1,lwd=1+as.numeric(is_S))

s3 <- rExp(pp, return_val="eigenvector")
s3s <- sort(s3)
plot(seq_along(s3s),s3s,type="n",log="y")
text(seq_along(s3s),s3s,labels=names(s3s))
## get_evec() is what calculates the eigenvector
## get_evec needs a testify flag! or a ...

s4 <- sort(make_state(N=1e6,E0=1000, params=pp))
s4 <- s4[s4>0]
plot(seq_along(s4),s4,type="n",log="y")
text(seq_along(s4),s4,labels=names(s4))


## make_ratemat

state <- make_state(pp[["N"]],E0=pp[["E0"]],params=pp)
M <- make_ratemat(state,pp)

## brute-force beta quantile (inverse CDF) function
## params as in qbeta()
Qbeta0 <- function(p,shape1,shape2,lower.tail=FALSE) {
  fn <- function(x) {pbeta(x,shape1,shape2,lower.tail=lower.tail)-p}
  uniroot(fn,interval=c(0,1))$root
}
## vectorized version 
Qbeta <- Vectorize(Qbeta0,c("p","shape1","shape2"))

## qbeta() works sometimes, Qbeta() works sometimes ...
## try qbeta(), switch to Qbeta() if it throws a warning
## vectorized manually because we need to detect warnings at
## the individual qbeta() evaluation level
Qbeta2 <- function(p,shape1,shape2,lower.tail=FALSE) {
  n <- max(length(p),length(shape1),length(shape2))
  p <- rep(p,length.out=n)
  shape1 <- rep(shape1,length.out=n)
  shape2 <- rep(shape2,length.out=n)
  res <- rep(NA,n)
  for (i in seq(n)) {
    res[i] <- tryCatch(qbeta(p[i],shape1[i],shape2[i],lower.tail=lower.tail),
                       warning=function(w) {
                         Qbeta0(p[i],shape1[i],shape2[i],lower.tail=lower.tail)
                       })
  }
  res
}

##' @param i prevalence
##' @param t testing proportion
##' @param phi dispersion
##' @param numint use numerical integration?
##' @param qfun beta quantile function to use
##' @param plot.it illustrate?
##' @param phiscale ??
## FIXME: remove vestigial (qfun, numint) options
## FIXME: test phi=0; add special case if necessary (and maybe phi=1?)
prop_pos_test0 <- function(i,t,phi,
                           numint=FALSE,qfun=Qbeta2,
                           plot.it=FALSE,
                           phiscale="constrained",
                           debug=FALSE
                           ) {
    if (debug) cat("p1: ",i,t,phi,"\n")  
    ## cat(phi,"\n")
    if (phiscale=="constrained") {
        ## transform from (0,inf) to (0,1)
        ## y = inverse(1-exp(-x)) -> -log(1-y)
       phi_0 <- phi
       phi <- -log(1-phi)
       if (is.nan(phi)) {
        cat(phi_0)
        }
    }
  
    i<- min(i,1) # treat i>1 as i=1 (exp growth of I makes i>1? or is it just the optimizer doing its thing?)
    
    a <- i/phi; b <- (1-i)/phi
    if (plot.it) curve(dbeta(x,a,b),from=0,to=1)
    ## need special logic for extreme phi values?
    if (debug) cat("p2: ",t,a,(1-i)/phi,"\n")
    lwr <- qfun(t,a,b,lower.tail=FALSE)
    if (debug) cat("lwr: ",lwr,"\n")
    if (plot.it) abline(v=lwr,col=2,lty=2)
    ## replace with pbeta(...,lower.tail=FALSE) with appropriate multiplier?
    val <- if (numint) {
               fn <- function(x) x*dbeta(x,a,b)
               integrate(fn, lower=lwr, upper=1)$value
           } else {
               pbeta(lwr,a+1,b,lower.tail=FALSE)*(a/(a+b))
           }
    if (debug) cat("val: ",lwr,"\n")
    return(val/pbeta(lwr,a,b,lower.tail=FALSE))
}
prop_pos_test <- Vectorize(prop_pos_test0,c("i","t","phi"))


## 'pure' version of prop_pos_test0, for RTMB use
prop_pos_test1 <- function(i,t,phi) {
    ## transform phi??
    a <- i/phi; b <- (1-i)/phi
    lwr <- qbeta(t,a,b,lower.tail=FALSE)
    val <- pbeta(lwr,a+1,b,lower.tail=FALSE)*(a/(a+b))
    return(val/pbeta(lwr,a,b,lower.tail=FALSE))
}


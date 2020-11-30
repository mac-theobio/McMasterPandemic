## What are good curve shapes for bridging?

m <- seq(0, 12, length.out=200)

## Saturating (attains start point)
S <- 2
F <- 0.5
h <- 6
sat <- (F-S)*m/(h+m) + S
plot(m, sat)

## Exponential (attains start point)
C <- 6
x <- exp(-m/C)
ex <- x*S + (1-x)*F
plot(m, ex)

## Logistic (attains nothing)
z <- 6
X <- 3
l <- plogis((m-z)/X)
lo <- (1-l)*S + l*F
plot(m, lo)


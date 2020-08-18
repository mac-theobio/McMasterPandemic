library(McMasterPandemic)
library(testthat)

## devtools::load_all("../..")
pp <- read_params("PHAC_testify.csv")
##  state <- make_state(N=1, E0=1e-5, type="ICU1")  ## FIXME: don't assume ICU1?
r1 <- rExp(params=pp, return_val="sim")
## including state breaks this ... ?
matplot(r1[,1],r1[,2:10],log="y",type="l",lty=1)
r2 <- rExp(params=pp, return_val="sim",testify=TRUE)
## plot S only
## matplot(r2[,1],r2[,2:5],log="y",type="l",lty=1)
## matplot(r2[,1],r2[,6:41],log="y",type="l",lty=1)
r2[nrow(r2),6:41]

library(tidyverse)
r2L <- (r2
    %>% as_tibble()
    %>% rename(date="t")
    %>% select(-c(D,X,N,P,foi))
    %>% McMasterPandemic:::pivot.pansim()
    %>% separate(var, into=c("pref","test"), sep="_")
)
gg1 <- ggplot(r2L,aes(date,value))+geom_line(aes(linetype=test))+facet_wrap(~pref) +
    scale_y_log10()
print(gg1)

print(gg1 %+% filter(r2L,date<30))


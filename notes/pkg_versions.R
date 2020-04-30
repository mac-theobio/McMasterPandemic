library(tidyverse)
get_pkgs <- function(fn="test0_ML.Rout") {
    ss <- scan(file=fn,what=character(1))
    pkgs <- grep("^[[:alpha:]]+_[0-9.-]+$",ss,value=TRUE)
    ss2 <- strsplit(pkgs,"_")
    return(tibble(pkg=vapply(ss2,"[[",character(1),1),
                  ver=vapply(ss2,"[[",character(1),2))
           %>% arrange(pkg)
           )
}

all_pkgs <- purrr::map(list(BB="test0.Rout",JD="test0_JD.Rout",ML="test0_ML.Rout"), get_pkgs)
comb_pkgs <- (Reduce(function(x,y) merge(x,y,by="pkg",all=TRUE),
                     all_pkgs)
    %>% as_tibble()
    %>% setNames(c("pkg",c("BB","JD","ML")))
    ## %>% filter(JD!=BB | ML!=BB)
    %>% filter(ML!=BB)
)

library(tidyverse)
library(ggplot2)
library(pins)

## MIDAS estimates retrieved from
midas_url <- "https://raw.githubusercontent.com/midas-network/COVID-19/master/parameter_estimates/2019_novel_coronavirus/estimates.csv"
x_temp <- (read_csv(pin(midas_url)))
x2<-(x_temp %>% select(authors,name,location_name,value_type,
              value,uncertainty_type,lower_bound, upper_bound,
              population)
## fix typos/inconsistencies
%>% mutate(name=tolower(name),
           name=gsub("hospitilization","hospitalization",name),
           name=gsub("reproductive","reproduction",name))
%>% right_join(focal)
## force numeric (FIXME, track down NAs later)
%>% mutate_at(c("value","lower_bound","upper_bound"), as.numeric)
%>% mutate(authora=collapse_authors(authors))
)

write.csv(as.data.frame(x2),"../inst/params/midas_estimates_sep.csv", row.names = FALSE)


focal <- read.csv(header=TRUE, text="
shortname, name
R0,basic reproduction number
inv_gamma_m1,time from symptom onset to isolation
inv_gamma_m2,time from symptom onset to recovery
inv_gamma_s1,time from symptom onset to hospitalization
inv_sigma,incubation period
latent,latent period
serial,serial interval
generation,generation interval
beta0,transmission rate")

## lastname of first author + initials of others (could be improved;
##  perhaps just disambiguating first author (Qifang+1, Qifang+2)
##  would be sufficient
collapse_authors <- function(a) {
    s <- strsplit(a,",")
    sapply(s,
           function(x) paste(c(gsub("^([[:alpha:]]+) .*","\\1",x[1]),
                               "+",
                               sapply(x[-1],function(y) substr(trimws(y),1,1))),
                              collapse=""))
}

## assumes we are in notes/ subdirectory ...
x <- (read_csv("../inst/params/midas_estimates.csv")
    %>% select(authors,name,location_name,value_type,
               value,uncertainty_type,lower_bound, upper_bound,
               population)
    ## fix typos/inconsistencies
    %>% mutate(name=tolower(name),
               name=gsub("hospitilization","hospitalization",name),
               name=gsub("reproductive","reproduction",name))
    %>% right_join(focal)
    ## force numeric (FIXME, track down NAs later)
    %>% mutate_at(c("value","lower_bound","upper_bound"), as.numeric)
    %>% mutate(authora=collapse_authors(authors))
)

    
gg0 <- (ggplot(x,aes(
              y=substr(authora,1,10), ## still not short enough
              x=value,xmin=lower_bound,xmax=upper_bound))
    + facet_wrap(~name,scale="free")
)

gg0 + geom_pointrange()

## when there are multiple values for an author, what do they represent?

## what about missing bounds?
## what about negative lower bounds?
## among-sample variance **not** incorporated yet.
##   (mean(var_comb) + var_among ??? this would be more appropriate for
##    multiple measurements of the _same_ data, which seems more appropriate here)
##   ( var(inv-var-wtd) mean without uncertainty is harmonic mean of var)

combine_values <- function(mid, lwr, upr) {
    logsd <- log(upr/lwr)/(2*qnorm(0.975))
    tau <- 1/logsd^2
    logmu <- log(mid)
    logmu_comb <- sum(logmu*tau)/sum(tau)
    ## ???
    ## logsd_comb <- sqrt(1/sum(tau)) ## doesn't use among-sample var!
    logsd_comb <- sqrt(mean(logsd^2) + var(logmu))
    data.frame(logmu=logmu_comb,logsd=logsd_comb)
}

comb <- (x  %>% drop_na(value,lower_bound,upper_bound)
    %>% filter(lower_bound>0)
    %>% group_by(name)
    %>% do(combine_values(.$value,.$lower_bound,.$upper_bound))
    %>% mutate(value=exp(logmu),
               lower_bound=exp(logmu-1.96*logsd),
               upper_bound=exp(logmu+1.96*logsd))
    %>% mutate(authora="zzzcombined")
)

xx_comb <- bind_rows(dplyr::select(comb,-c(logmu,logsd)),
                     dplyr::select(x,
                                   authora,value,
                                   lower_bound,upper_bound,name))
gg0 %+%  xx_comb + geom_pointrange(aes(colour=(authora=="zzzcombined"))) +
    scale_colour_manual(values=c("black","red"),guide=FALSE)

ggsave("midas_estimates.pdf")

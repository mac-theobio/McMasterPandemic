## fitting log-linear models to death/hospital/ICU time series
library(tidyverse)
library(ggstance)
library(broom)
library(broom.mixed)
library(ggplot2); theme_set(theme_bw())
library(colorspace)
library(glmmTMB)
library(lme4)

x <- readRDS("data/20200402-location_covid_data.rds")
## Location, Date, Admitted, Hospitalized, ICU, Deaths, Discharged

x2 <- (x
    %>% pivot_longer(names_to="var",-c(Location,Date))
    %>% na.omit()
)
gg1 <- (ggplot(x2,aes(Date,value,colour=Location))
    + geom_line()
    + geom_point()
    + facet_wrap(~var,scale="free")
    + theme(legend.position="none")
)
## linear, all vars
print(gg1)
vars <- c("Hospitalized","ICU","Deaths")
## interesting vars only
print(gg1F <- gg1 %+% filter(x2,var %in% vars))
## log scale
print(gg1F + scale_y_log10())

## could do quadratic fit and predict slope at end of time period??
## deaths definitely steeper across the board

fitfun <- function(var="Hospitalized",family="nbinom2") {
    cat(var,"\n")
    x3 <- (x2[x2$var==var,] ## ugh, non-tidy
        %>% mutate(t=as.numeric(Date),
                   t0=t-min(t))
    )
    cat(format(min(x$Date)),"\n")
    nb2 <- glmmTMB(value~t0 + (t0|Location),
                   family=family,
                   data=x3)
    return(nb2)
}

## NB fits are problematic for ICU and Deaths, so use Poisson for now
family <- c("nbinom2","poisson","poisson")
fitList <- purrr::map2(vars,family,fitfun) %>% setNames(vars)
                      
## distribution of slopes
nb2 <- fitList$Hospitalized
hist(ranef(nb2)$cond$Location[,"t0"],main="",
     xlab="Deviation from population growth rate (/day)")
## caterpillar plot??

## nb2_mer <- glmer.nb(value~t0 + (t0|Location),
##                     data=x3)
## ## intercept is useful but we'd like to be able to see slope alone!
## lattice::dotplot(ranef(nb2_mer))

## nb1 <- update(nb2, family="nbinom1")
## poiss <- update(nb2, family="poisson")

## (coef(nb2)$cond$Location
##     %>% rownames_to_column("Location")
##     %>% as_tibble() %>% arrange(desc(t0))
## )

tidy1 <- function(fit) {
    res <- (tidy(fit,"ran_vals")
        %>% filter(term=="t0")
        %>% select(-c(effect,component,group,term))
        %>% mutate(lwr=estimate-std.error,upr=estimate+std.error)
        %>% rename(location="level")
    )
    return(res)
}


tt1 <- map_dfr(fitList,tidy1,.id="var")
c_order <- (filter(tt1, var=="Hospitalized")
    %>% arrange(estimate)
    %>% pull(location)
)
tt1 <- mutate(tt1, location=factor(location,levels=c_order))

print(ggplot(tt1,aes(y=location,x=estimate,xmin=lwr,xmax=upr))
      + geom_pointrange()
      + facet_wrap(~var)
      )
## FIXME: figure out missing/NA values?

tidy2 <- function(fit) {
    res <- (tidy(fit,"ran_vals")
        %>% dplyr::select(-c(effect,component,group))
        %>% mutate(lwr=estimate-std.error,upr=estimate+std.error)
        %>% dplyr::select(-std.error)
        %>% rename(location="level")
        %>% mutate_at("term",~ifelse(.=="(Intercept)","int",.))
        %>% pivot_wider(names_from="term",
                        values_from=c("estimate","lwr","upr"))
    )
    return(res)
}

tt2 <- map_dfr(fitList,tidy2,.id="var")

## y = log(2)/x 
trans_dbltime <- scales::trans_new("dbltime",
                                   function(x) log(2)/x, function(x) log(2)/x,
                                   domain=c(0,Inf))

## saturation or controls?
print(ggplot(tt2,aes(estimate_int,
                     estimate_t0,colour=var,shape=var))
      + geom_point(size=2)
      + geom_linerange(aes(xmin=lwr_int,xmax=upr_int))
      + geom_linerange(aes(ymin=lwr_t0,ymax=upr_t0))
      + scale_colour_discrete_qualitative()
      + labs(y="deviation of growth rate (/day)",
             x="deviation of log(predicted number on 25 Mar)")
      )

## FIXME: better scales, coeffs instead of ranefs


## ugh: abbreviate names
names(fitList) <- ifelse(names(fitList)=="ICU","ICU",substr(names(fitList),1,1))
cc <- (map_dfr(fitList,~tibble::rownames_to_column(coef(.)$cond$Location),
               .id="var")
    %>% setNames(c("var","location","int","slope"))
)
saveRDS(cc,file="data/location_slopes.rds")

library(EpiEstim)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
source("notes/ontario_clean.R")
incid <- (ont_recent
    %>% filter(var=="newConfirmations")
    %>% pull(value)
)


## see ontario_calibration.Rmd
schoolClose <- "2020-03-17"
countryClose <- "2020-03-23"
socialClose <- "2020-03-28"
bd <- as.Date(c(schoolClose,countryClose,socialClose))
e1 <- estimate_R(incid=incid,
           method="uncertain_si",
           config=make_config(list(mean_si=6,
                                   std_mean_si = 1,
                                   min_mean_si = 3,
                                   max_mean_si = 8,
                                   std_si=1.5,
                                   std_std_si=0.5,
                                   min_std_si=0.5,
                                   max_std_si=2.5,
                                   n1=100,n2=100)))
plot(e1)

res <- (e1$R
    %>% rename_at(vars(starts_with("Q")),
                  ~gsub("Quantile\\.0\\.([0-9]+)\\(R\\)","q\\1",.))
    %>% rename(med="Median(R)")
    ## FIXME: not 100% sure how dates line up?
    %>% mutate(date=head(sort(unique(ont_recent$date)),nrow(e1$R)))
)
mle_cal <- (readRDS("git_push/ontario_R0_cal1.rds")
    %>% rename(date="start_date")
    ## HACK: scale_x_date has no oob argument?
    %>% mutate(date=c(min(res$date)-5,as.Date(date[-1])))
    %>% bind_rows(tibble(date=max(ont_recent$date),
                         r0=tail(.$r0,1),
                         R0=tail(.$R0,1),
                         Gbar=tail(.$Gbar,1),
                         dbl_time=tail(.$dbl_time,1)))
)

print(gg1 <- ggplot(res,aes(date, y=med))
      + geom_line()
      + geom_ribbon(aes(ymin=q05,ymax=q95),colour=NA,alpha=0.5,fill="red")
      + geom_ribbon(aes(ymin=q025,ymax=q975),colour=NA,alpha=0.2,fill="blue")
      + scale_y_log10(breaks=c(0.75, 1, 1.1, 1.3, 1.5, 2, 3, 5))
      + geom_vline(xintercept=bd,lty=2)
      + geom_step(data=mle_cal,aes(y=R0),
                  lwd=2,alpha=0.3)
      )


## q25/q75 from EpiEstim are nonsense, not sure why?
## CIs on 


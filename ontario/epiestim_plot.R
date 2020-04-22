library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

f_args <- ont_cal1$forecast_args
ss <- summary(ont_cal1)
mle_cal0 <- (ss
    %>% rename(date="start_date")
    ## HACK: scale_x_date has no oob argument?
    ## set first date to a more recent time
    %>% mutate(date=c(min(epiestim_fit$date)-5,as.Date(date[-1])))
    %>% bind_rows(tibble(date=max(ont_recent$date),
                         r0=tail(.$r0,1),
                         R0=tail(.$R0,1),
                         Gbar=tail(.$Gbar,1),
                         dbl_time=tail(.$dbl_time,1)))
)
    
mle_cal1 <- (mle_cal0
    %>% full_join(tibble(date=seq.Date(min(.$date),max(.$date),by="1 day")),by="date")
    %>% arrange(date)
    %>% select(date,R0)   ## other columns too?
    %>% tidyr::fill(R0)
    ## OK to use base_params (we haven't calibrated anything about reporting delay)
    %>% mutate(R0conv=calc_conv(R0,update(f_args$base_params,c_prop=1)))
)

gg1 <- (ggplot(epiestim_fit,aes(date, y=med))
      + geom_line()
      + geom_ribbon(aes(ymin=q05,ymax=q95),colour=NA,alpha=0.5,fill="red")
      + geom_ribbon(aes(ymin=q025,ymax=q975),colour=NA,alpha=0.2,fill="blue")
      + scale_y_log10(breaks=c(0.75, 1, 1.1, 1.3, 1.5, 2, 3, 5))
      + geom_vline(xintercept=bd,lty=2)
      + geom_line(data=mle_cal1,aes(y=R0conv),
                  lwd=2,alpha=0.3)

)
print(gg1)

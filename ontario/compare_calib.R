library(McMasterPandemic)
library(tidyverse)
library(directlabels)
library(bbmle)
library(ggplot2); theme_set(theme_bw())
dl_offset <- 25

## files loaded via wrapR

cal_list <- nlist(ont_cal1,ont_cal_1brks,ont_cal_2brks,hosponly=ont_cal2,
                  ont_cal_noICU, ont_cal_noICU_1brks, ont_cal_noICU_2brks,
                  ont_cal_HD_2brks,
                  ont_cal_noICU_2brks_prior)

## abbreviate names
names(cal_list) <- gsub("ont_cal_?","",names(cal_list))
names(cal_list)[names(cal_list)=="1"] <- "3brks"

keep_vars <- c("death","H","ICU","report")
pp <- (map_dfr(cal_list,predict,.id="cal")
    %>% dplyr::filter(var %in% keep_vars
             , date>as.Date("2020-03-01"))
)
gg1 <- (ggplot(pp, aes(date,value,colour=var))
    + geom_line(aes(lty=cal))
    + facet_wrap(~var,scales="free")
    + scale_y_log10(limits=c(1,NA))
    + geom_point(data=filter(trans_state_vars(ont_all), var %in% keep_vars))
    + geom_dl(method="last.bumpup",aes(label=cal))
    + expand_limits(x=max(pp$date)+dl_offset)
)
ggsave(plot=gg1,"compare_calib.1.pdf",width=10,height=7)

pp2 <- (pp
    %>% filter(grepl("noICU|HD",cal))
    %>% mutate_at("cal",~gsub("noICU_","",.))
)
gg2 <- gg1 %+% pp2
ggsave(plot=gg2,"compare_calib.noICU.pdf",width=10,height=7)

fit <- ont_cal_noICU_2brks_prior
summary(fit)
sqrt(diag(vcov(fit$mle2)))
## pop_pred_samp(fit$mle2,PDify=TRUE,return_wts=TRUE,n=200)

## stuff below will fail if vcov is not PD ...
try({

    pred1 <- predict(fit,ensemble=TRUE)
    plot(pred1)

    pred1S <- predict(fit,ensemble=TRUE,
                   stoch=c(proc=TRUE,obs=TRUE),
                   new_params=list(obs_disp=20,proc_disp=1),
                   end_date="1-July-2020") 
 print(plot(pred1S)
    + geom_line(data=pred1,lty=2)
     + geom_point(data=filter(trans_state_vars(ont_all), var %in% keep_vars))
 )
})


library(McMasterPandemic)
library(tidyverse)
library(directlabels)
library(bbmle)
load("ontario_calibration.RData")  ## ont_cal1, 3 breaks, all types
load("ontario_calibration_1brks.RData")  ## ont_cal_1brks
load("ontario_calibration_2brks.RData")  ## ont_cal_2brks
load("ontario_calibration_hosponly.RData") 
load("ontario_calibration_noICU.RData") 
load("ontario_calibration_noICU_1brks.RData")
load("ontario_calibration_noICU_2brks.RData")
load("ontario_calibration_hosponly.RData")
load("ontario_calibration_HD_2brks.RData")
load("ontario_clean.RData")

cal_list <- nlist(ont_cal1,ont_cal_1brks,ont_cal_2brks,hosponly=ont_cal2,
                  ont_cal_noICU, ont_cal_noICU_1brks, ont_cal_noICU_2brks,
                  ont_cal_HD_2brks)
keep_vars <- c("death","H","ICU","report")
pp <- (map_dfr(cal_list,predict,.id="cal")
    %>% dplyr::filter(var %in% keep_vars
             , date>as.Date("2020-03-01"))
)
(gg1 <- ggplot(pp, aes(date,value,colour=var))
    + geom_line(aes(lty=cal))
    + facet_wrap(~var,scales="free")
    + scale_y_log10(limits=c(1,NA))
    + geom_point(data=filter(trans_state_vars(ont_all), var %in% keep_vars))
    + geom_dl(method="last.bumpup",aes(label=cal))
    + expand_limits(x=as.Date("2020-05-01"))
)

print(gg1 %+% filter(pp,grepl("noICU|HD",cal)))

fit <- ont_cal_noICU_2brks
summary(fit)
sqrt(diag(vcov(fit$mle2)))
pred1 <- predict(fit,ensemble=TRUE)
plot(pred1)

pred1S <- predict(fit,ensemble=TRUE,
                  stoch=c(proc=TRUE,obs=TRUE),
                  new_params=list(obs_disp=20,proc_disp=1),
                  end_date="1-July-2020") 
(plot(pred1S)
    + geom_line(data=pred1,lty=2)
    + geom_point(data=filter(trans_state_vars(ont_all), var %in% keep_vars))
)

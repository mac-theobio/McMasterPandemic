library(McMasterPandemic)
library(tidyverse)
theme_set(theme_bw())
library(ggforce)
## load("ont_cal_splinecomp.RData")
## L <- load(".ont_keep.RData")

## should have saved this ..
inputs <- expand.grid(spline_df=3:7,
                      spline_sd_pen=c(NA,0.001,0.01,0.1),
                      knot_quantile_var=c(NA,"H","report"))

nm <- apply(inputs,1,paste,collapse="_")
names(res_list) <- nm
all_pred <- (purrr::map_dfr(res_list, predict, .id="mod")
    %>% separate(mod,into=c("df","pen","svar"),sep="_")
    %>% filter(var %in% c("H","death","report"))
    %>% mutate_at(c("pen","svar"), ~ ifelse(.=="NA" | is.na(.), "none", .))
)
var_fun <- function(rvar,dat=ont_all_sub, pred_dat=all_pred, date_lims=as.Date(c(NA,NA))) {
    if (is.na(date_lims[1])) date_lims[1] <- min(pred_dat$date)
    if (is.na(date_lims[2])) date_lims[2] <- max(pred_dat$date)
    fdate <- . %>% filter(var==rvar, between(date, date_lims[1], date_lims[2]))
    pred_df <- pred_dat %>% fdate()
    dat <- full_join(dat %>% fdate(),
                     pred_df %>% select(svar,df) %>% distinct(),
                     by=character())
    gg1 <- (ggplot(pred_df,aes(date,value))
        + facet_grid(svar~df, labeller=label_both)
        + geom_line(aes(colour=factor(pen)))
        + scale_y_log10()
        + scale_x_date(limits=date_lims)
        + ggtitle(sprintf("var=%s",rvar))
        + geom_point(data=dat,alpha=0.2)
        + theme(panel.spacing=grid::unit(0,"lines"))
    )
    return(gg1)
}

var_fun("H")
var_fun("H", date_lims=as.Date(c("2020-04-15",NA)))
var_fun("report")
var_fun("report", date_lims=as.Date(c("2020-04-15",NA)))

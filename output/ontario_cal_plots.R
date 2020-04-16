library(McMasterPandemic)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(directlabels)

keep_vars <- c("H","ICU","d","incidence","report","newTests/1000")
f_args <- attr(ont_cal1,"forecast_args")

get_type <- . %>%  mutate(vtype=ifelse(var %in% c("incidence","report","d"),
                                       "inc","prev"))
sub_vars <- . %>% dplyr::filter(var %in% keep_vars)

## make forecast with specified end date
mk_fc <- function(new_end=NULL) {
    if (!is.null(new_end)) {
        f_args$end_date <- new_end
    }
    fc <- (do.call(forecast_sim,
                      c(list(p=ont_cal1$par), f_args))
        %>% sub_vars()
        %>% get_type()
    )
    return(fc)
}

## pivot and unpivot to rescale newTests ...
ont_all_sub <- (ont_all
    %>% mutate_at("var",trans_state_vars)
    %>% pivot_wider(names_from="var",values_from="value",id_cols="date")
    %>% mutate(`newTests/1000`=newTests/1000)
    %>% select(-newTests)
    %>% pivot_longer(names_to="var",-date)
    %>% get_type()
)
## data frame for labeling new tests
newtest_lab <-data.frame(date=as.Date("2020-04-10"),
                         value=10,
                         vtype="prev",
                         var="newTests/1000",
                         lab="newTests/1000")
## data frame for labeling ICU capacities
capac_info <- data.frame(value=c(630,1300),
                         vtype="prev",
                         var="ICU",
                         lab=c("current","expanded"))


plotfun <- function(forecast, data=ont_all_sub, breaks=bd,
                    dlspace=1,
                    limspace=10,
                    add_tests=TRUE,
                    add_ICU_cap=FALSE) {
    data <- sub_vars(data)
    if (!add_tests) data <- filter(data,var!="newTests/1000")
    p <- (ggplot(forecast,aes(date,value,colour=var))
        + geom_line()
        + scale_y_log10(limits=c(1,NA),oob=scales::squish)
        + geom_point(data=data)
        + geom_line(data=data,alpha=0.2)
        + geom_vline(xintercept=breaks,lty=2)
        + facet_wrap(~vtype,ncol=1,scale="free_y")
        + labs(y="")
        + theme(legend.position="none",
                ## https://stackoverflow.com/questions/10547487/remove-facet-wrap-labels-completely
                strip.background = element_blank(),
                strip.text = element_blank())
    )
    ## trying to fix spacing on the fly: kluge!
    ## dlspace not found
    ## geom_dl() must do weird evaluation env stuff
    if (add_tests) {
        p <- p + geom_text(data=newtest_lab,aes(label=lab))
    }
    if (add_ICU_cap) {
        p <- (p
            + geom_hline(data=capac_info,aes(yintercept=value,
                                             colour=var),lty=3)
            + geom_text(data=capac_info,aes(y=value,x=min(data$date),
                                            label=lab),vjust=-1)
        )
    }
    return(p + geom_dl(method=list(dl.trans(x=x+1),cex=1,'last.bumpup'),
                       aes(label=var))
           + expand_limits(x=max(forecast$date)+limspace)
           )
}

forecast1 <- mk_fc()
gg_cal1 <- plotfun(forecast1,limspace=15)
ggsave(gg_cal1, file="output/ont_cal1.pdf")
ggsave(gg_cal1, file="output/ont_cal1.png")


forecast2 <- mk_fc(new_end="2020-05-15")
gg_cal2 <- plotfun(forecast2,add_tests=FALSE,add_ICU_cap=TRUE,limspace=20)
ggsave(gg_cal2, file="output/ont_cal2.pdf")
ggsave(gg_cal2, file="output/ont_cal2.png")

forecast3 <- mk_fc(new_end="2020-10-01")
gg_cal3 <- plotfun(forecast3,add_tests=FALSE,limspace=45)
ggsave(gg_cal3, file="output/ont_cal3.pdf")
ggsave(gg_cal3, file="output/ont_cal3.png")


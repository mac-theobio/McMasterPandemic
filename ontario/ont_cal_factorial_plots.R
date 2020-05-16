library(McMasterPandemic)
library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(colorspace)
L <- load("ont_cal_factorial.RData")

ss <- function(flags,i) as.logical(as.numeric(substr(flags,i,i)))
get_term <- function(flags,coefs) {
    terms <- character(0)
    if (ss(flags,2))  terms <- c(terms,sprintf("mobility (power=%1.2f)",
                                               coefs$time_beta[1]))
    if (ss(flags,3)) terms <- c(terms,sprintf("spline (bs, %d knots)",7))
    if (ss(flags,4)) terms <- c(terms,sprintf("phenom_het (power=%1.2f)",
                                              coefs$params[["zeta"]]))
    if (length(terms)==0) terms <- "nothing"
    return(paste(terms,collapse="\n"))
}
predfun <- function(r) {
    predict(r$fit) %>% filter(var %in% c("H","death","report"))
    ## FIXME: restrict to dates for each var? right_join etc.
}

bad <- sapply(res_list, inherits, "try-error")
terms <- map2(factorial_combos[!bad],res_list[!bad],
              ~get_term(.x, coef(.y$fit,"fitted")))
factorial_combos[bad]  ## phenom_het only; mobility + phenom_het
res_list <- setNames(res_list[!bad], terms)

pp <- purrr::map_dfr(res_list, predfun, .id="model")
gg1 <- (ggplot(pp,aes(date,value,colour=var,shape=var,lty=var))
    + geom_line()
    + facet_wrap(~model)
    + scale_y_log10()
    ## data is identical for all facets
    + geom_point(data=r_list[[1]]$data,alpha=0.5)
    ## limit dates to those with available data
    + scale_x_date(limits=range(r_list[[1]]$data$date))
    + scale_colour_discrete_qualitative()
)
ggsave(gg1,file="fac_plot.pdf")

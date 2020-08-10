library(ggplot2);theme_set(theme_bw()+theme(panel.spacing=grid::unit(0,"lines")))
library(dplyr)
library(tidyr)
library(cowplot)

source("makestuff/makeRfuns.R")
commandEnvironments()
if (!interactive()) {
    makeGraphics()
} else {
    load("testify_sim.rda")
}

ymin <- 1    
ymax <- 1e6  ## screen out pathological values

## ugh: modify test positivity to percentage so visible
simdat <- (simdat
    %>% pivot_wider(names_from=var,values_from=value)
    %>% mutate_at("positivity", ~ . * 100)
    %>% pivot_longer(-c(date, W_asymp, iso_t, testing_intensity), names_to="var")
)

totaltest_data <- (tibble(
    var="total_test",
    testing_intensity=unique(simdat$testing_intensity))
    %>%     mutate(value=testing_intensity*params[["N"]])
)
    
gg <- (ggplot(simdat)
    + aes(x=date,y=value)
    + geom_line()
    + scale_color_manual(values=c("black","red","blue","orange","purple"))
    + scale_y_log10(limits=c(ymin, ymax),oob=scales::squish)
    + ylab("Daily count")
    + theme(legend.position = "bottom")
)


ff <- function(i,data=simdat) filter(data, testing_intensity==i)
mm <- function(i) ggtitle(sprintf("testing intensity=%1.2g",i))
for (i in unique(simdat$testing_intensity)) {
    ggx <- (gg
        %+% ff(i)
        + mm(i)
        + geom_hline(data=ff(i,totaltest_data),aes(yintercept=value),lty=2)
    )
    print(ggx + facet_grid(W_asymp~iso_t, labeller=label_both) + aes(color=var))
    print(ggx + facet_grid(W_asymp~var, labeller=label_both) + aes(color=factor(iso_t)))
    print(ggx + facet_grid(iso_t~var, labeller=label_both) + aes(color=factor(W_asymp)))
}

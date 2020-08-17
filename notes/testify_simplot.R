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

simdat <- (simdat
    %>% pivot_wider(names_from=var,values_from=value)
    %>% pivot_longer(-c(date, testing_type, iso_t, omega), names_to="var")
)

totaltest_data <- (tibble(
    var="total_test",
    omega=unique(simdat$omega))
    %>%     mutate(value=omega*params[["N"]])
)
    
gg <- (ggplot(simdat)
    + aes(x=date,y=value)
    + geom_line()
    + scale_color_manual(values=c("black","red","blue","orange","purple"))
    + scale_y_log10(limits=c(ymin, ymax),oob=scales::squish)
    + ylab("Daily count")
    + theme(legend.position = "bottom")
)


ff <- function(i,data=simdat) filter(data, omega==i)
mm <- function(i) ggtitle(sprintf("omega=%1.2g",i))
for (i in unique(simdat$omega)) {
    ggx <- (gg
        %+% ff(i)
        + mm(i)
        + geom_hline(data=ff(i,totaltest_data),aes(yintercept=value),lty=2)
    )
    print(ggx + facet_grid(testing_type~iso_t, labeller=label_both) + aes(color=var))
    print(ggx + facet_grid(testing_type~var, labeller=label_both) + aes(color=factor(iso_t)))
    print(ggx + facet_grid(iso_t~var, labeller=label_both) + aes(color=factor(testing_type)))
}

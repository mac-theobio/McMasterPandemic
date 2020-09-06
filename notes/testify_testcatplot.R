library(ggplot2);theme_set(theme_bw()+theme(panel.spacing=grid::unit(0,"lines")))
library(dplyr)
library(tidyr)
library(cowplot)
library(colorspace)
scale_colour_discrete <- colorspace::scale_colour_discrete_qualitative

source("makestuff/makeRfuns.R")
commandEnvironments()
if (!interactive()) {
    makeGraphics()
} else {  ## this is out dated?
    load("testify_sim.rda")
}

ymin <- 1    
ymax <- 1e6  ## screen out pathological values

print(simdat)


gg <- (ggplot(simdat)
    + aes(x=date,y=value)
    + geom_line(aes(linetype=testcat))
    + scale_y_log10()
)

ff <- function(i,data=simdat) filter(data, testing_intensity==i)
mm <- function(i) ggtitle(sprintf("testing intensity=%1.2g",i))

## FIXME: now we have two more dimensions to deal with (pref or testcat), that's annoying
for (i in unique(simdat$testing_intensity)) {
    ggx <- (gg
        %+% ff(i)
        + mm(i)
    )
    print(ggx + facet_grid(W_asymp~iso_t, labeller=label_both) + aes(color=pref))
    print(ggx + facet_grid(W_asymp~pref, labeller=label_both) + aes(color=factor(iso_t)))
    print(ggx + facet_grid(iso_t~pref, labeller=label_both) + aes(color=factor(W_asymp)))
}

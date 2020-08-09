library(ggplot2);theme_set(theme_bw()+theme(panel.spacing=grid::unit(0,"lines")))
library(dplyr)
library(tidyr)
library(cowplot)

source("makestuff/makeRfuns.R")
commandFiles()

## ugh: modify test positivity to percentage so visible
simdat <- (simdat
    %>% pivot_wider(names_from=var,values_from=value)
    %>% mutate_at("positivity", ~ . * 100)
    %>% pivot_longer(-c(date, W_asymp, iso_t), names_to="var")
)

gg <- (ggplot(simdat)
	+ aes(x=date,y=value)
	+ geom_line()
	+ scale_color_manual(values=c("black","red","blue","orange","purple"))
	+ scale_y_log10(limits=c(0.1, NA))
	+ ylab("Daily count")
	+ theme(legend.position = "bottom")
)

print(gg + facet_grid(W_asymp~iso_t, labeller=label_both) + aes(color=var))
print(gg + facet_grid(W_asymp~var, labeller=label_both) + aes(color=factor(iso_t)))
print(gg + facet_grid(iso_t~var, labeller=label_both) + aes(color=factor(W_asymp)))

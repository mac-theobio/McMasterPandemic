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

scientific_10 <- function(x,suppress_ones=TRUE) {
   s <- scales::scientific_format()(x)
   ## substitute for exact zeros
   s[s=="0e+00"] <- "0"
   ## regex: [+]?  = "zero or one occurrences of '+'"
   s2 <- gsub("e[+]?", " %*% 10^", s )
   ## suppress 1 x
   if (suppress_ones) s2 <- gsub("1 %\\*% +","",s2)
   parse(text=s2)
}

varnames <- data.frame(var = c("total_test", "report", "postest", "pos_per_million", "incidence")
	, varname = c("Daily Test", "Report", "Positive Test", "Positive test per million", "Incidence")
)


ymin <- 1    
ymax <- 1e6  ## screen out pathological values

print(simdat)


simdat <- (simdat
	 %>% left_join(.,varnames)
	 %>% mutate(Isolation = ifelse(iso_t == 1, "Yes","No")
	 	, Gbar = paste0("Gbar"," = ",Gbar," days")
		, Gbar = factor(Gbar, levels=c("Gbar = 6 days", "Gbar = 12 days"))
		)
	 %>% filter(var != "report")
	 %>% select(-c(iso_t,testing_type,omega,var))
	 %>% pivot_wider(names_from=varname,values_from=value)
    %>% pivot_longer(-c(date, Gbar, Isolation, testing_intensity), names_to="var")
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
    + scale_y_log10(labels = scientific_10, limits=c(ymin, ymax),oob=scales::squish)
    + ylab("Daily count")
    + theme(legend.position = "bottom")
)


ff <- function(i,data=simdat) filter(data, testing_intensity==i)
mm <- function(i) ggtitle(sprintf("testing intensity=%1.2g",i))
for (i in unique(simdat$testing_intensity)[1]) {
    ggx <- (gg
#       %+% ff(i)
#		  %+% mm(i)
#       + geom_hline(data=ff(i,totaltest_data),aes(yintercept=value),lty=2)
    )
#    print(ggx + facet_grid(Gbar~iso_t, labeller=label_both) + aes(color=var))
    print(ggx + facet_grid(Gbar~var) + aes(color=Isolation))
#    print(ggx + facet_grid(iso_t~var, labeller=label_both) + aes(color=factor(Gbar)))
}

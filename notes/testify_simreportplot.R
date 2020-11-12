library(ggplot2);theme_set(theme_bw()+theme(panel.spacing=grid::unit(0,"lines")))
library(dplyr)
library(tidyr)
library(cowplot)

callArgs <- "testwt_N.simplots.Rout testify_simplot.R testwt_N.sims.rda"
source("makestuff/makeRfuns.R")
commandEnvironments()
if (!interactive()) makeGraphics()

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

varnames <- data.frame(var = c("total_test", "positivity", "postest", "CumIncidence", "incidence")
	, varname = c("Daily Test", "Positivity", "Positive tests\nper day", "Cumulative Incidence", "Incidence")
)


simdat2 <- (simdat
	 %>% filter(iso_t %in% c(0,1))
	 %>% filter(omega == 0.25)
	 %>% filter(testing_intensity == 0.01)
    %>% left_join(.,varnames)
    %>% mutate(Isolation = ifelse(iso_t == 1, "Yes","No")
             , speed = factor(Gbar, levels=c(6,12), labels=c("faster","slower"))
               )
    %>% filter(var != "report")
    %>% select(-c(iso_t,testing_type,omega,var))
    %>% pivot_wider(names_from=varname,values_from=value)
#    %>% mutate(`% positive tests`= `Positive test per million`/1e4)
    %>% select(-c(`Daily Test`, Gbar))
    %>% pivot_longer(-c(date, speed, Isolation, testing_intensity), names_to="var")
)


ymin <- 1    
ymax <- NA

print(simdat2)

gg <- (ggplot(simdat2)
    + aes(x=date,y=value,colour=Isolation, size=Isolation)
    + geom_line(alpha=0.7)
    + scale_color_manual(values=c("black","red"))
    + scale_size_manual(values=c(2,1))
    ## + scale_y_log10(labels = scientific_10, limits=c(ymin,ymax), oob=scales::squish)
    + ylab("Daily count")
    + theme(legend.position = "bottom")
    + facet_grid(speed~var)
)

gg1 <- gg %+% filter(simdat2, var=="% positive tests") + labs(y="Percent")
gg2 <- gg %+% filter(simdat2, var!="% positive tests") + scale_y_log10(limits=c(1,NA))

# plot_grid(gg2,gg1,nrow=1,rel_widths=c(2,1))



ggall <- (ggplot(simdat)
			 + aes(x=date,y=value,colour=factor(iso_t), linetype=factor(omega))
			 + scale_colour_manual(values=c("red","blue","black"))
			 + geom_line()
			 + facet_grid(var~testing_intensity, scale="free")
)

print(ggall 
	%+% (simdat 
		%>% filter(var != "total_test")
		%>% filter(Gbar == 6)
		)
	+ ggtitle(targetname())
)


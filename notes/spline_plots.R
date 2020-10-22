library(ggplot2);theme_set(theme_bw())
library(dplyr)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(combo_pars)

print(spline_df)

combo_pars <- (combo_pars
	%>% mutate(spline_pen = ifelse(grepl("pen",mod),"yes","no")
		, spline_pen = ifelse(type == "true", "no", spline_pen)
		, model = ifelse(grepl("E0",mod),"withE0","withoutE0")
	)
)


gg <- (ggplot(combo_pars, aes(x=mod, y=value,color=mod))
	+ geom_point()
	+ facet_wrap(~var,scale="free")
	+ scale_alpha_manual(values=c(1,0.3,0.3,0.3,0.3))
	+ scale_color_manual(values=c("blue","blue","red","red","black"))
)

print(gg)

print(spline_df)

ggsplines <- (ggplot(spline_df,aes(time,color=mod,group=seed))
	+ facet_wrap(~mod)
	+ scale_color_manual(values=c("blue","blue","red","red","black"))
)

print(ggsplines + geom_line(aes(y=bt)))
print(ggsplines + geom_line(aes(y=Rt)))


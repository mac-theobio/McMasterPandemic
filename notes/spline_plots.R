library(ggplot2);theme_set(theme_bw())

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(combo_pars)

gg <- (ggplot(combo_pars, aes(x=mod, y=value,color=mod, alpha=type))
	+ geom_point()
	+ facet_wrap(~var,scale="free")
	+ scale_alpha_manual(values=c(1,0.3,0.3))
	+ scale_color_manual(values=c("black","blue","red"))
)

print(gg)

print(spline_df)

ggsplines <- (ggplot(spline_df,aes(time,bt,color=mod,alpha=mod,group=seed))
	+ geom_line()
	+ facet_wrap(~mod)
	+ scale_color_manual(values=c("black","blue","red"))
	+ scale_alpha_manual(values=c(1,0.3,0.3))
)

print(ggsplines)

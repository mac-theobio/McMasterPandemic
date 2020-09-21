library(ggplot2);theme_set(theme_bw())

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(combo_pars)

gg <- (ggplot(combo_pars, aes(x=var, y=value,color=type, alpha=type))
	+ geom_point()
	+ facet_wrap(~var,scale="free")
	+ scale_alpha_manual(values=c(0.3,1))
	+ scale_color_manual(values=c("red","black"))
)

print(gg)

print(spline_df)

ggsplines <- (ggplot(spline_df,aes(time,bt,color=type,alpha=type,group=seed))
	+ geom_line()
	+ scale_color_manual(values=c("red","black"))
	+ scale_alpha_manual(values=c(0.3,1))
)

print(ggsplines)

library(ggplot2);theme_set(theme_bw())

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(combo_pars)

print(spline_df)

gg <- (ggplot(combo_pars, aes(x=mod, y=value,color=spline_pen, alpha=type, shape = spline_pen))
	+ geom_point()
	+ facet_wrap(~var,scale="free")
	+ scale_alpha_manual(values=c(1,0.3,0.3,0.3,0.3))
	+ scale_color_manual(values=c("black","orange","purple"))
)

print(gg)

print(spline_df)

ggsplines <- (ggplot(spline_df,aes(time,color=mod,alpha=mod,group=seed))
	+ facet_grid(spline_pen~mod)
	+ scale_color_manual(values=c("black","blue","red","orange","purple"))
	+ scale_alpha_manual(values=c(1,0.3,0.3,0.3,0.3))
)

print(ggsplines + geom_line(aes(y=bt)))
print(ggsplines + geom_line(aes(y=Rt)))


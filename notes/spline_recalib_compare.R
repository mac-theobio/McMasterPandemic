library(McMasterPandemic)
library(tidyverse);theme_set(theme_bw())

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(plot(ff, data=dat))

aa <- coef(ff, "fitted")
bb <- coef(ff_refit, "fitted")
spline_df <- data.frame(time = 1:nrow(X)
	, true_spline = exp(X[,-1] %*% matrix(aa$time_beta, ncol=1))
	, refit_spline = exp(X[,-1] %*% matrix(bb$time_beta, ncol=1))
)

print(spline_df)

sdf <- (spline_df
	%>% gather(key="type",value="value",-time)
)

gg <- (ggplot(sdf,aes(x=time,y=value,color=type))
	+ geom_line()
	+ scale_color_manual(values=c("red","black"))
)

print(gg)


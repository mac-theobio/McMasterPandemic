library(ggplot2);theme_set(theme_bw())
library(dplyr)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

combo_pars$withb0 <- "yes"
combo_ss_pars$withb0 <- "no"

combo_pars <- bind_rows(combo_pars, combo_ss_pars)

true_pars <- combo_pars %>% filter(type == "true")
combo_pars <- combo_pars %>% filter(type != "true")

spline_df2 <- bind_rows(spline_df, spline_shape_df)

print(spline_df)

combo_pars <- (combo_pars
	%>% mutate(spline_pen = ifelse(grepl("pen",mod),"yes","no")
		, spline_pen = ifelse(type == "true", "no", spline_pen)
		, model = ifelse(grepl("E0",mod),"withE0","withoutE0")
		, model = ifelse(grepl("spline_shape",seed),paste0(model,"withoutB0"),paste0(model,"withB0"))
		, model = ifelse(model == "withoutE0withoutB0", "none"
			, ifelse(model == "withoutE0withB0", "B0"
				, "E0B0")
			)
	)
)

spline_df2 <- (spline_df2
	%>% mutate(model = ifelse(grepl("spline_shape",seed),"noB0","B0"))
)

gg <- (ggplot(combo_pars, aes(x=mod, y=value,color=model))
	+ geom_point(position= position_dodge(width = 0.9))
	+ facet_wrap(~var,scale="free")
	+ geom_hline(data=true_pars, aes(yintercept=value))
	+ scale_alpha_manual(values=c(1,0.3,0.3,0.3,0.3))
	+ scale_color_manual(values=c("black","red","blue"))
	+ coord_flip()
)

print(gg)

print(unique(spline_df$mod))
print(unique(spline_df2$mod))


ggsplines <- (ggplot(spline_df2,aes(time,color=mod,group=seed))
	+ facet_grid(model~mod)
#	+ scale_color_manual(values=c("blue","red"))
	+ ylim(c(0,10))
)

# print(ggsplines + geom_line(aes(y=bt)))
print(ggsplines + geom_line(aes(y=Rt)))


library(ggplot2);theme_set(theme_bw())
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

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




print(head(dat))

print(tail(dat))

varnames <- data.frame(var=c("hosp","death","positivity","postest","total_test")
	, varname = c("Hospital Admission", "Death", "Positivity","Positive Test", "Daily Test")
)

dat <- left_join(dat,varnames)

gg <- (ggplot(dat, aes(x=date))
	+ geom_line(aes(y=value))
	+ geom_point(aes(y=data),alpha=0.4)
	+ facet_wrap(~varname, scale="free", ncol=2)
#	+ scale_y_log10(limits=c(1,NA))
	+ scale_y_log10(labels = scientific_10)
	+ xlim(c(as.Date("2020-04-01"),as.Date("2020-08-01")))
	+ ylab("Daily Counts")
	+ xlab("Date")
	+ ggtitle("Logistic testing")
)

print(gg)



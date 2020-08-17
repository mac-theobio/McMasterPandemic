library(tidyverse);theme_set(theme_bw())
library(McMasterPandemic)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

## bad. repeating code 

print(comboframe)

ln <- list.files(pattern = "RDS", path = "./cachestuff")

print(ln)

clean_res <- function(x){

	res <- readRDS(paste0("cachestuff/",x))
	if(class(res$fit) == "NULL"){return(data.frame())}
	input <- unlist(strsplit(x,split="[.]"))[2]

	dd <- predict(res$fit
		, ensemble=FALSE
		, keep_vars=c("postest","death","H")
	)
	ddcombo <- (res$data
		%>% transmute(date, postest, death, H)
		%>% gather(key = "var", value="data",-date)
		%>% left_join(dd,.)
		%>% mutate(testing_intensity = comboframe[input,1]
			, keep_vars = comboframe[input,"keep_vars"]
			, opt_testify = comboframe[input,"opt_testify"]
			, testing_type = comboframe[input,"testing_type"]
		)
	)

	return(ddcombo)
}

res_list <- lapply(ln,clean_res)

dat <- bind_rows(res_list)

saveVars(dat)

library(tidyverse);theme_set(theme_bw())
library(McMasterPandemic)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

## bad. repeating code 

testing_intensity <- c(0.001, 0.01, 0.1)
keep_vars <- c("postest", "H/death", "postest/H/death")
opt_testify <- c(TRUE,FALSE)

comboframe <- expand.grid(testing_intensity=testing_intensity
        , keep_vars = keep_vars
        , opt_testify = opt_testify
)



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
			, keep_vars = comboframe[input,2]
			, opt_testify = comboframe[input,3]
		)
	)

	return(ddcombo)
}

res_list <- lapply(ln,clean_res)

dat <- bind_rows(res_list)

saveVars(dat)

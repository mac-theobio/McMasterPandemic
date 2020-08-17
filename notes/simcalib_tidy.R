library(tidyverse);theme_set(theme_bw())
library(McMasterPandemic)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

## bad. repeating code 

comboframe <- expand.grid(keep_vars = keep_vars
   , opt_testify = opt_testify
   , testing_type = testing_type
)

print(comboframe)

# ln <- list.files(pattern = "RDS", path = "./cachestuff")

# print(ln)

ln <- "simcalib.1.RDS"

clean_res <- function(x){

	res <- readRDS(paste0("cachestuff/",x))
	if(class(res$fit) == "NULL"){return(data.frame())}
	input <- unlist(strsplit(x,split="[.]"))[2]

	dd <- predict(res$fit
		, ensemble=FALSE
		, keep_vars=c("postest","negtest","hosp","death")
	)
	ddcombo <- (res$data
		%>% transmute(date, postest, death, hosp, negtest)
		%>% gather(key = "var", value="data",-date)
		%>% left_join(dd,.)
		%>% mutate(testing_intensity = comboframe[input,1]
			, keep_vars = comboframe[input,"keep_vars"]
			, opt_testify = comboframe[input,"opt_testify"]
			, testing_type = comboframe[input,"testing_type"]
		)
	)

	## stupid way to back calculate total_test
	postest <- (ddcombo 
		%>% filter(var == "postest")
	)
	negtest <- (ddcombo
		%>% filter(var == "negtest")
	)

	total_test <- data.frame(date = postest$date
		, var = "total_test"
		, value = postest$value + negtest$value
		, data = postest$data + negtest$data
	)
	positivity <- data.frame(date = postest$date
		, var = "positivity"
		, value = postest$value/total_test$value
		, data = postest$data/total_test$data
	)
	ddcombo2 <- (ddcombo
		%>% filter(var != "negtest")
		%>% bind_rows(.,total_test)
		%>% bind_rows(.,positivity)
	)

	return(ddcombo2)
}

res_list <- lapply(ln,clean_res)

dat <- bind_rows(res_list)

saveVars(dat)

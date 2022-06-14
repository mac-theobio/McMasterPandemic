library(ggplot2)
library(McMasterPandemic)
library(tidyr)

params = read_params("PHAC_testify.csv")

testified_model = make_testified_model(params=params, start_date="2020-01-01", end_date="2021-01-01")
testified_history = simulation_history(testified_model)

(testified_history
  %>% select(Date, any_of(names(testified_model$state)))
  %>% pivot_longer(-Date)
  %>% separate(name, into = c("epi", "test"), sep = "_",remove = FALSE)
  %>% ggplot()
  + facet_grid(epi~., scales = "free_y")
  + geom_line(aes(Date, value, colour=test))
)


























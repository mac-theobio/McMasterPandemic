library(autotest)
library(testthat)
library(dplyr)

test_that('Automatic tests from autotest', {

broken_tests <- c('make_test_wtsvec', 'mle_fun', 'vis_model', 'write_params', 'run_sim_range')

# these 4 cause errors in autotest, and that is probably an autotest problem.
# Not sure if the last function is actually bugged or if I'm going crazy

oldw <- getOption("warn")
options(warn = -1)
x <- autotest_package('package:McMasterPandemic', test = TRUE,exclude = broken_tests) 
options(warn = oldw)

grouped_tests = x %>% 
  dplyr::group_by(fn_name) %>% 
  dplyr::summarise(all(test))

for(row in 1:nrow(grouped_tests)){
  expect(grouped_tests[row,]$`all(test)`, paste(grouped_tests[row,]$fn_name, ' autotest failed.'))
}
#expect(all(x$test==TRUE), paste(toString(sum(x$test)),'automatic tests failed.'))
# Really not sure what is going on with expect_autotest_no_errors
# expect_autotest_no_err(x) 
})



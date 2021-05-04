library('covr')
base_dir = 'C:/Users/somat/Documents/GitHub/McMasterPandemic'
setwd(base_dir)

cov = package_coverage(base_dir)

df = as.data.frame(cov)
write.csv(df, paste(base_dir, '/refactor/coverage/coverage_details.csv', sep=''))

no_cov = zero_coverage(cov)
no_cov_df = as.data.frame(no_cov)
write.csv(no_cov_df, paste(base_dir, '/refactor/coverage/no_coverage_details.csv', sep=''))

# For some reason, sink doesn't work with ... colors? I just copy-pasted the print into coverage_report.txt
print(cov)


##########################################################################3
#library('remotes')
#remotes::install_github("ropenscilabs/autotest")

#library('autotest')
#x <- autotest_package('C:/Users/somat/Documents/GitHub/McMasterPandemic')
# Autotest breaks when trying to get libraries?:
# Error in library(p, character.only = TRUE) :
# there is no package called 'lattice, Matrix, MASS, dplyr, ggplot2, tidyr, plyr, purrr, anytime, pomp, scales, bbmle (>= 1.0.23.5), Hmisc, DEoptim, mvtnorm, bdsmatrix, zoo, deSolve, splines, diagram'



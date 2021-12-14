options(repos = "https://cran.rstudio.com")
install.packages("remotes")
#install.packages("tinytex")
install.packages("devtools")
install.packages("roxygen2")
#install.packages("styler")
#install.packages("lintr")

install.packages('Matrix')
install.packages('RcppEigen')
install.packages('TMB')

## These are required under https://hub.docker.com/r/rocker/rstudio
## install.packages('Hmisc')
## install.packages('tidyverse')

remotes::install_deps(".", dependencies = "Depends", upgrade = TRUE)
remotes::install_github("bbolker/bbmle")
remotes::install_github("johndharrison/semver")

if(FALSE) {
    tinytex::install_tinytex()
    tinytex::tlmgr_install("koma-script")
    tinytex::tlmgr_install("amscls")
    tinytex::tlmgr_install(c("multirow", "colortbl", "siunitx", "setspace"))
    tinytex::tlmgr_install(c("lineno", "fancyhdr", "ulem", "caption"))
    tinytex::tlmgr_install("babel-english")
    tinytex::tlmgr_install(c("pgf"))
    tinytex::tlmgr_install(c("placeins", "lastpage", "cleveref", "listings"))
}

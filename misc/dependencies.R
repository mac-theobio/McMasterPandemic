install.packages('remotes')
remotes::install_deps('.', dependencies = TRUE, upgrade = TRUE)

remotes::install_github('bbolker/bbmle')

install.packages('tinytex')
tinytex::install_tinytex()
tinytex::tlmgr_install('koma-script')
tinytex::tlmgr_install('amscls')
tinytex::tlmgr_install(c('multirow','colortbl','siunitx','setspace'))
tinytex::tlmgr_install(c('lineno','fancyhdr','ulem','caption'))
tinytex::tlmgr_install('babel-english')
tinytex::tlmgr_install(c('pgf', 'preview', 'xcolor'))
tinytex::tlmgr_install(c('placeins','lastpage','cleveref','listings'))

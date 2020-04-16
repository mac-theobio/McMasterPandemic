## This is McMaster pandemic: ontario subdir

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

vim_session:
	bash -cl "vmt"

######################################################################

Sources += $(wildcard */*.R)

Sources += $(wildcard man/*.Rd) NAMESPACE

######################################################################

ontario_nbfit.Rout: ontario_clean.Rout ontario_nbfit.R
	$(run-R)

ontario_calibration.Rout:  ontario_clean.Rout ontario_calibration.R
ontario_cal_plots.Rout: ontario_clean.Rout ontario_calibration.Rout ontario_cal_plots.R

epiestim_plot.Rout: ontario_clean.Rout ontario_calibration.Rout epiestim.Rout epiestim_plot.R

notes/ontario_calibration_report.html: ontario_clean.Rout ontario_calibration.Rout epiestim.RData epiestim_plot.Rout ontario_cal_plots.Rout ontario_nbfit.Rout.png ont_cal1.png notes/ontario_calibration_report.Rmd

%.html: %.Rmd
	Rscript -e 'library("rmarkdown"); render("$<", output_format="html_document")'


### Makestuff

Sources += Makefile

## Sources += content.mk
## include content.mk

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls $@

localstuff: 
	ln -s ../makestuff .
	ls $@

-include makestuff/os.mk

-include makestuff/wrapR.mk
-include makestuff/rmd.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk

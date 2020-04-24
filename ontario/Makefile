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

pkgdesc = $(../DESCRIPTION)

######################################################################

ontario_nbfit.Rout: ontario_clean.Rout ontario_nbfit.R

ontario_calibration.Rout:  ontario_clean.Rout ontario_calibration.R $(pkgdesc)

ontario_calibration_hosponly.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_hosponly.R

ontario_calibration_2brks.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_2brks.R

ontario_calibration_1brks.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_1brks.R

ontario_calibration_noICU.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_noICU.R

ontario_calibration_noICU_1brks.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_noICU_1brks.R

ontario_calibration_noICU_2brks.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_noICU_2brks.R

ontario_calibration_HD_2brks.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_HD_2brks.R

ontario_calibration_2brks_ndt.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_2brks_ndt.R

ontario_calibration_noICU_2brks_prior.Rout:  ontario_clean.Rout ontario_calibration.Rout ontario_calibration_noICU_2brks_prior.R


compare_calib.Rout:  ontario_clean.Rout ontario_calibration.RData ontario_calibration_noICU_1brks.RData ontario_calibration_noICU_2brks.RData ontario_calibration_noICU.RData ontario_calibration_2brks.RData ontario_calibration_1brks.RData ontario_calibration_hosponly.RData ontario_calibration_HD_2brks.RData ontario_calibration_noICU_2brks_prior.RData

epiestim.Rout: ontario_clean.Rout
epiestim_plot.Rout: ontario_clean.Rout ontario_calibration.Rout epiestim.Rout epiestim_plot.R

ontario_calibration_report.html: ontario_clean.RData ontario_calibration.RData epiestim.RData epiestim_plot.Rout.png ontario_nbfit.Rout.png ont_cal1.png ontario_calibration_report.Rmd

McMasterReport.pdf: McMasterReport.Rnw ontario_calibration.RData ontario_calibration_2brks.RData ontario_calibration_hosponly.RData


%.html: %.Rmd
	Rscript -e 'rmarkdown::render("$<", output_format="html_document")'

%.pdf: %.Rnw
	Rscript -e 'knitr::knit2pdf("$<")'

### Makestuff

Sources += Makefile

## Sources += content.mk
## include content.mk

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	ln -s ../makestuff .
	ls $@

-include makestuff/os.mk

-include makestuff/wrapR.mk
-include makestuff/rmd.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk

# McMasterPandemic

<!-- badges: start -->
[![R-CMD-check](https://github.com/mac-theobio/McMasterPandemic/workflows/R-CMD-check/badge.svg)](https://github.com/mac-theobio/McMasterPandemic/actions)
<!-- badges: end -->

Compartmental epidemic models for forecasting and analysis of infectious disease pandemics: contributions from Ben Bolker, Jonathan Dushoff, David Earn, Morgan Kain, Michael Li, Irena Papst (in alphabetical order). Feedback is welcome at the [issues list](https://github.com/mac-theobio/McMasterPandemic/issues), or e-mail us.

### Documentation
* See the accompanying [pkgdown](https://mac-theobio.github.io/McMasterPandemic/) site.
* Functions are under the `Reference` tab.
* Vignettes (under construction!) are under `Articles`.
* [Issues](https://github.com/mac-theobio/McMasterPandemic/issues) list.
* [Shiny app](https://mcmasterpandemic.shinyapps.io/mcmasterpandemicshiny/) written by Zach Levine

### Installation

The repository contains an R package and various workflows/analyses. You can fork/clone the repository (from [here](https://github.com/mac-theobio/McMasterPandemic)) and install locally or use `remotes::install_github("mac-theobio/McMasterPandemic")` to install the package. You will need to first install the developer version of `bbmle` (`remotes::install_github("bbolker/bbmle")`) before installing `McMasterPandemic`.

### For developers

* to re-install the package, including re-building and incorporating vignettes, use `make build`
* If you modify function arguments, you should change the roxygen documentation accordingly. If you change the roxygen documentation, please use `make doc-update` to update the `.Rd` files.
* **please test/check the package periodically** as you go (use `make pkgcheck` and `make pkgtest` from the shell or `devtools::check()` and `devtools::test()` from within R). (Tests are also run on [GitHub Actions](https://github.com/mac-theobio/McMasterPandemic/actions); if you want to skip CI testing, e.g. for a trivial commit, put `[skip ci]` somewhere in your commit message.) Please don't make a habit of pushing without testing.
* To avoid whitespace diffs between versions/branches, automatically style the package with `make style` or run `misc/macpan_style.R`. `make style` (or running `misc/macpan_lint.R`) also creates a new file, `misc/lints.csv`, which contains stylistic and other lints that `styler` cannot automatically fix.
* rebuild the `pkgdown` site using [GitHub Actions](https://github.com/mac-theobio/McMasterPandemic/actions/workflows/pkgdown.yaml): click "run workflow" on the link to rebuild. 
* Code that is used in the refactoring process should go in the top-level `refactor` folder. 
* Slow tests are skipped by default; this process is controlled by the presence of a `MACPAN_TEST_LEVEL` environment (shell) variable. Many of the test files contain this code:
```r
testLevel <- if (nzchar(s <- Sys.getenv("MACPAN_TEST_LEVEL"))) as.numeric(s) else 1
```
which sets `testLevel` to a default value of 1 unless the environment variable is found. You can set this environment variable outside of your R session (e.g. via `export MACPAN_TEST_LEVEL=2` in  a bash shell), or via `Sys.setenv(MACPAN_TEST_LEVEL=2)` from within R. In principle this mechanism allows for a hierarchy of slowness; at present only 1 and >1 are distinguished.

### Documentation 

The documentation is a little bit scattered right now, working on cleaning it up. In addition to the standard short descriptions of the functions (`help(package="McMasterPandemic")`), stuff can be found: 

* in the vignettes (look at the source code in the [vignettes] directory or `vignette(<title>, package="McMasterPandemic")`)
    * `getting_started`
	* `model`: design decisions and information for developers
	* `calibration` (**very out of date**)
	* `farr`: stuff on Farr's law and phenomenological curve-fitting (**very incomplete and likely to remain so for now**)
	* `testing_flow`: incorporating testing dynamics (**ditto**)
* `McMasterReport.Rnw`: this is a more or less up-to-date description of calibration to Ontario data
* `ontario_calibration_report.html`: more technical and less up-to-date than the preceding document
* `TODO.md`: active to-do list

More bits and pieces: `notes/refactor.Rmd`, `testing.md`, `reimplementation.md`

### DISCLAIMER

All use of this package is at your own risk. Quantitative forecasts are only as good as their parameter estimates.


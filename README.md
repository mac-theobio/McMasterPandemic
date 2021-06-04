# McMasterPandemic

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbolker/McMasterPandemic/workflows/R-CMD-check/badge.svg)](https://github.com/bbolker/McMasterPandemic/actions)
[![DOI](https://zenodo.org/badge/252492971.svg)](https://zenodo.org/badge/latestdoi/252492971)
<!-- badges: end -->


Compartmental epidemic models for forecasting and analysis of infectious disease pandemics. Feedback is welcome at the [issues list](https://github.com/bbolker/McMasterPandemic/issues).

### Installation

The repository contains an R package and various workflows/analyses. You can fork/clone the repository (from [here](https://github.com/bbolker/McMasterPandemic)) and install locally or use `remotes::install_github("bbolker/McMasterPandemic")` to install the package. You will need to first install the developer version of `bbmle` (`remotes::install_github("bbolker/bbmle")`) before installing `McMasterPandemic`. 

### Documentation 

See the accompanying [pkgdown site](https://bbolker.github.io/McMasterPandemic).

- Functions are under the `Reference` tab.
- Vignettes (under construction!) are under `Articles`
- [Shiny app](https://mcmasterpandemic.shinyapps.io/mcmasterpandemicshiny/) written by Zach Levine

### For developers

* to re-install the package, including re-building and incorporating vignettes, use `make build`
* re-build the `pkgdown.extras` site by tagging a commit message with `rebuild site`. For example, `Update README.md [skip ci][rebuild site]`.
* Automatically style the package with `make style` or run `misc/macpan_style.R`. Additionally, `make style` (or running `misc/macpan_lint.R` creates a new file, `misc/lints.csv`, which contains stylistic and other lints that styler cannot automatically fix.
* it's OK to push *small* changes to the master branch; use your own branch + pull request for anything non-trivial
* If you modify function arguments, please update the roxygen documentation and use `make doc-update` (or `devtools::document`) to update the `.Rd` files.
* **Check the package periodically** as you go (use `make pkgcheck` and `make pkgtest` from the shell or `devtools::check()` and `devtools::test()` from within R). (Tests are also run on [GitHub Actions](https://github.com/bbolker/McMasterPandemic/actions); if you want to skip CI testing, e.g. for a trivial commit, put `[skip ci]` somewhere in your commit message.) Don't push without checking.
* Code that is used in the refactoring process should go in the top-level `refactor` folder. 
* Slow tests are skipped by default; this process is controlled by the presence of a `MACPAN_TEST_LEVEL` environment (shell) variable. Many of the test files contain this code:
```r
testLevel <- if (nzchar(s <- Sys.getenv("MACPAN_TEST_LEVEL"))) as.numeric(s) else 1
```
which sets `testLevel` to a default value of 1 unless the environment variable is found. You can set this environment variable outside of your R session (e.g. via `export MACPAN_TEST_LEVEL=2` in  a bash shell), or via `Sys.setenv(MACPAN_TEST_LEVEL=2)` from within R. In principle this mechanism allows for a hierarchy of slowness; at present only 1 and >1 are distinguished.

## Obsolete documentation

(Keeping links here for now)

* in the vignettes (look at the source code in the [vignettes] directory or `vignette(<title>, package="McMasterPandemic")`)
    * `getting_started`
	* `model`: design decisions and information for developers
	* `calibration` (**very out of date**)
	* `testing_flow`: incorporating testing dynamics (**ditto**)
* `McMasterReport.Rnw`: this is a more or less up-to-date description of calibration to Ontario data
* `ontario_calibration_report.html`: more technical and less up-to-date than the preceding document
* `TODO.md`: active to-do list
* issues list on github

More bits and pieces: `notes/refactor.Rmd`, `testing.md`, `reimplementation.md`;shiny app [GitHub repo](https://github.com/ZachLevine-11/McMasterPandemicShiny).



### DISCLAIMER

All use of this package is at your own risk. Quantitative forecasts are only as good as their parameter estimates.


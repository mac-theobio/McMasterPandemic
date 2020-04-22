# McMasterPandemic

Compartmental epidemic models for forecasting and analysis of infectious disease pandemics: contributions from Ben Bolker, Jonathan Dushoff, David Earn, Morgan Kain, Michael Li (in alphabetical order). Feedback is welcome at the [issues list](https://github.com/bbolker/McMasterPandemic/issues), or e-mail us.

### Installation

The repository contains an R package and various workflows/analyses. You can fork/clone the repository (from [here](https://github.com/bbolker/McMasterPandemic)) and install locally or use `remotes:install_github("bbolker/McMasterPandemic")` to install the package.

### For developers

* to re-install the package, including re-building and incorporating vignettes, use `make build`
* If you modify function arguments, you should change the roxygen documentation accordingly. If you change the roxygen documentation, please use `make doc-update` to update the `.Rd` files.
* please test/check the package periodically as you go (use `make pkgcheck` and `make pkgtest` from the shell or `devtools::check()` and `devtools::test()` from within R).

### DISCLAIMER

All use of this package is at your own risk. Quantitative forecasts are only as good as their parameter estimates.


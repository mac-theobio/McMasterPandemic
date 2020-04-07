# McMasterPandemic

Compartmental epidemic models for forecasting and analysis of infectious disease pandemics.

Please send any comments to
- Ben Bolker <bolker@mcmaster.ca>
- David Earn <earn@math.mcmaster.ca>

### Installation

`remotes:install_github("bbolker/McMasterPandemic")` (you have to install the `remotes` package first, if you don't have it)

### For developers

* to re-install the package, including re-building and incorporating vignettes, use `make build`
* If you modify function arguments, you should change the roxygen documentation accordingly. If you change the roxygen documentation, please use `make doc-update` to update the `.Rd` files.
* please test/check the package periodically as you go (use `make pkgcheck` and `make pkgtest` from the shell or `devtools::check()` and `devtools::test()` from within R).

### DISCLAIMER

We expect this package to be valuable primarily for _qualitative_
understanding, which can help to guide planning during a pandemic.
Quantitative forecasts are only as good as the parameter estimates,
which are very challenging to make.


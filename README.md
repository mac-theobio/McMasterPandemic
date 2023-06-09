# McMasterPandemic

<hr>

> **Note**
> The McMasterPandemic project is [evolving](https://canmod.net/misc/macpan2_presentation)
> * All new development work is happening in the [macpan2 repository](https://github.com/canmod/macpan2).    
> * The [McMasterPandemic repository](https://github.com/mac-theobio/McMasterPandemic) is in **maintenance mode** to support **current users**. 
> * If you are a **new user**, thank you for your interest! Please start [here at macpan2](https://github.com/canmod/macpan2).

<hr>

<!-- badges: start -->
[![R-CMD-check](https://github.com/mac-theobio/McMasterPandemic/workflows/R-CMD-check/badge.svg)](https://github.com/mac-theobio/McMasterPandemic/actions)
[![Task list badge](https://img.shields.io/static/v1.svg?label=kanban&message=tmb%20engine&color=blue)](https://github.com/mac-theobio/McMasterPandemic/projects/7)

<!-- badges: end -->



Compartmental epidemic models for forecasting and analysis of infectious disease pandemics: contributions from Ben Bolker, Jonathan Dushoff, David Earn, Weiguang Guan, Morgan Kain, Michael Li, Irena Papst, Steve Walker (in alphabetical order). Feedback is welcome at the [issues list](https://github.com/mac-theobio/McMasterPandemic/issues), or e-mail us.

### Refactoring for Speed and Generality

We are [currently refactoring](https://github.com/mac-theobio/McMasterPandemic/projects/7) McMasterPandemic to make it faster and more general. We have merged the [development branch for this refactoring project](https://github.com/mac-theobio/McMasterPandemic/tree/tmb-condense) back into the master branch. However the classic functionality of McMasterPandemic is still available and coexists with the more general interface and faster engine. To get started with the classic functionality, please read [this vignette](https://mac-theobio.github.io/McMasterPandemic/articles/getting_started.html). To get started with the faster and more general functionality, please read [this user guide](https://canmod.github.io/macpan-book/).

### Documentation
* [Package documentation](https://mac-theobio.github.io/McMasterPandemic/)
* [User guide](https://canmod.github.io/macpan-book/)
* [Issues list](https://github.com/mac-theobio/McMasterPandemic/issues)
* [Roadmap](https://github.com/mac-theobio/McMasterPandemic/blob/master/TODO.md)
* [Draft contributors guide for developers](https://github.com/mac-theobio/McMasterPandemic/blob/master/CONTRIBUTING.md)

### Installation

The repository contains an R package and various workflows/analyses. This repository is not on [CRAN](https://cran.r-project.org/) so you will need to either fork/clone the repository (from [here](https://github.com/mac-theobio/McMasterPandemic)) or install directly from GitHub. Either option will (may?) require you to first install two R packages that are also not on CRAN.
```
remotes::install_github("bbolker/bbmle")
remotes::install_github("johndharrison/semver")
```
Note that these commands depend on having the `remotes` package, which you can get from CRAN with the following command from an R prompt.
```
install.packages('remotes')
```

To install `McMasterPandemic` itself you follow the same formula.
```
remotes::install_github("mac-theobio/McMasterPandemic")
```
If this command fails it may be because your R installation is not set up to compile C++ code. Windows users should be able to get past this issue by installing [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

The classic `McMasterPandemic` functionality described [here](https://mac-theobio.github.io/McMasterPandemic/articles/getting_started.html) does not depend on C++ code, and you can get access to this functionality by installing with this command.
```
remotes::install_github("mac-theobio/McMasterPandemic@v0.0.20.1")
```

Getting access to experimental features can sometimes be achieved with this command.

```
remotes::install_github("mac-theobio/McMasterPandemic@tmb-condense")
```



### DISCLAIMER

All use of this package is at your own risk. Quantitative forecasts are only as good as their parameter estimates.

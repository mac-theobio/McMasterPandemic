##
## vignette("rd",package="roxygen2")
##

## dealing with LaTeX in roxygen documentation:
## https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Conditional-text

##' \pkg{McMasterPandemic}
##'
##' This R package provides Compartmental epidemic models for
##' forecasting and analysis of infectious disease pandemics.
##' 
##' See the vignette \code{McMasterPandemic.pdf} for now...
##' Eventually \code{vignette("McMasterPandemic",package="McMasterPandemic")}
##' should work...
##'
##' While (even the strictly deterministic version of) the model is
##' not implemented as ODEs, its solutions are similar to those of the
##' following ODEs.
##' 
##' \ifelse{latex}{
##'   \out{
##'     \begin{array}[rlc]
##'       S' & = & - (\beta_0 / N) S (C_a I_a + C_p I_p + (1-iso_m)C_m I_m + (1-iso_s)C_s I_s) \\
##'       E' & = &  \\
##'       & \vdots &  \\
##'       R' & = &
##'     \end{array}
##'   }
##' }{
##'   \ifelse{html}{
##'     \out{
##'        d<i>S</i>/d<i>t</i> = &minus; (<i>&beta;</i><sub>0</sub>/<i>N</i>) <i>S</i> [<i>C</i><sub>a</sub> <i>I</i><sub>a</sub> + <i>C</i><sub>p</sub> <i>I</i><sub>p</sub> + (1&minus;iso<sub>m</sub>)<i>C</i><sub>m</sub> <i>I</i><sub>m</sub> + (1&minus;iso<sub>s</sub>)<i>C</i><sub>s</sub> <i>I</i><sub>s</sub>]
##'     }
##'   }{
##'   }
##' }
##'
##' \eqn{ dE/dt =  }
##' 
##' etc
##' 
##' \eqn{ dR/dt =  }
##'
##' Use \code{browseVignettes("McMasterPandemic")} to find vignettes
##' for this package.
"_PACKAGE"

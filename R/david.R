
##' Print method for \code{panpamams} objects
##'
##' @param x a \code{panparams} object
##'
##' @export
print.panparams <- function( x ) {
    param_meanings <- c(
        beta0 = "transmission rate",
        Ca = "relative asymptomatic transmissibility",
        Cp = "relative pre-symptomatic transmissibility",
        Cs = "relative severely symptomatic transmissibility",
        Cm = "relative mildly symptomatic transmissibility",
        alpha = "proportion of infections that are asymptomatic",
        gamma = "1 / mean LATENT period",
        lambda_a = "1 / mean days in asymptomatic infectious class",
        lambda_s = "1 / mean days in severely symptomatic infectious class",
        lambda_m = "1 / mean days in mildly symptomatic infectious class",
        lambda_p = "1 / mean days in pre-symptomatic infectious class",
        rho = "1 / mean days in acute care",
        delta = "proportion in acute care patients who die",
        mu = "proportion of symptomatic infections that are mild",
        N = "total population size",
        E0 = "number of initially exposed individuals",
        iso_m = "proportion of mildly symptomatic patients who are isolated",
        iso_s = "proportion of severely symptomatic patients who are isolated",
        phi1 = "proportion of severe infections that do NOT require ICU",
        phi2 = "proportion of ICU patients who die",
        psi1 = "1 / mean days in ICU if survive",
        psi2 = "1 / mean days in ICU if die",
        psi3 = "1 / mean days post-ICU until discharge"
    )
    x_meanings <- param_meanings[names(x)]
    xout <- data.frame(value=round(as.numeric(x),3),
                       meaning=x_meanings)
    knitr::kable(xout)
    return(invisiable(xout))
}

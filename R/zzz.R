.onLoad <- function(lib, pkg) {
    options(MP_badsum_action="stop",
            MP_badsum_tol=1e-12)
}

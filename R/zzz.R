.onLoad <- function(lib, pkg) {
    flex_version <- readLines(
        system.file("tmb/recommended_spec_version",
            package = "McMasterPandemic"
        )
    )[1]
    options(
        MP_badsum_action = "warning",
        MP_badsum_tol = 1e-12,
        MP_flex_spec_version = flex_version,
        MP_flex_spec_doc_site = "https://canmod.net/misc/flex_specs",
        MP_flex_spec_dll = "McMasterPandemic"
    )
}

.onUnload <- function(libpath) {
    library.dynam.unload("McMasterPandemic", libpath)
}

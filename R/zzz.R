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
        MP_flex_spec_dll = "McMasterPandemic",
        # https://stackoverflow.com/questions/8396577/check-if-character-value-is-a-valid-r-object-name
        MP_name_search_regex = "((([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*)|[.])",
        MP_use_state_rounding = TRUE
    )
}

.onUnload <- function(libpath) {
    library.dynam.unload("McMasterPandemic", libpath)
}

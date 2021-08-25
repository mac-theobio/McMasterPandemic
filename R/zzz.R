.onLoad <- function(lib, pkg) {
    options(
        MP_badsum_action = "warning",
        MP_badsum_tol = 1e-12,
        MP_flex_spec_version = "0.0.1",
        MP_flex_spec_doc_site = "https://canmod.net/misc/flex_specs"
    )
}
